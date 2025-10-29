#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import pysam
from Bio import SeqIO

def load_cb_umi_mapping(tsv_path):
    """
    Load a TSV with 4 columns:
        1) Raw CB
        2) Corrected CB
        3) Raw UMI
        4) Corrected UMI

    Returns:
        cb_map: dict[str raw_cb] = str corrected_cb
        umi_map_by_raw_cb: dict[str raw_cb] = dict[str raw_umi] = str corrected_umi
    """
    cb_map = {}
    umi_map_by_raw_cb = defaultdict(dict)

    open_func = gzip.open if tsv_path.endswith('.gz') else open
    mode = 'rt' if tsv_path.endswith('.gz') else 'r'

    line_ct = 0
    malformed = 0
    pairs_ct = 0
    unique_cb = set()

    try:
        with open_func(tsv_path, mode) as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 4:
                    malformed += 1
                    continue

                raw_cb, corr_cb, raw_umi, corr_umi = parts[0], parts[1], parts[2], parts[3]
                cb_map[raw_cb] = corr_cb
                umi_map_by_raw_cb[raw_cb][raw_umi] = corr_umi
                pairs_ct += 1
                unique_cb.add(raw_cb)
                line_ct += 1
    except Exception as e:
        print(f"Error loading CB/UMI mapping from {tsv_path}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded mapping from {tsv_path}")
    print(f"  Lines read: {line_ct + malformed}")
    print(f"  Malformed lines skipped (<4 cols): {malformed}")
    print(f"  Unique raw CBs: {len(unique_cb)}")
    print(f"  (raw CB, raw UMI) pairs: {pairs_ct}")
    return cb_map, umi_map_by_raw_cb


def parse_10x_barcode_from_r1(sequence, cb_length=16, umi_length=12):
    """
    Parse cell barcode and UMI from 10x R1 read sequence.
    Returns (cell_barcode, umi, remaining_sequence)
    """
    if len(sequence) < cb_length + umi_length:
        return None, None, ""
    cell_barcode = sequence[:cb_length]
    umi = sequence[cb_length:cb_length + umi_length]
    remaining_sequence = sequence[cb_length + umi_length:]
    return cell_barcode, umi, remaining_sequence


def process_pair(r1_record, r2_record, cb_map, umi_map_by_raw_cb, stats, cb_length=16, umi_length=12):
    """
    Process a single FASTQ pair and return (r1_bam, r2_bam), or None if skipped.
    Mapping logic:
      - Look up corrected CB by raw CB.
      - Look up corrected UMI by (raw CB, raw UMI).
    Both must exist to emit a pair.
    """
    # Normalize names (handle /1 and /2 suffixes if present)
    r1_name = r1_record.id.removesuffix('/1')
    r2_name = r2_record.id.removesuffix('/2')
    if r1_name != r2_name:
        stats['name_mismatch'] += 1
        return None

    r1_seq = str(r1_record.seq)
    cell_barcode, umi, remaining_r1_seq = parse_10x_barcode_from_r1(r1_seq, cb_length, umi_length)
    if cell_barcode is None or umi is None:
        stats['too_short'] += 1
        return None
    if len(remaining_r1_seq) == 0:
        stats['no_remaining_seq'] += 1
        return None

    # CB and UMI corrections (UMI in context of *raw* CB)
    corrected_cb = cb_map.get(cell_barcode)
    if corrected_cb is None:
        stats['cb_not_in_map'] += 1
        return None

    umi_map_for_cb = umi_map_by_raw_cb.get(cell_barcode)
    if not umi_map_for_cb:
        stats['cb_has_no_umi_map'] += 1
        return None

    corrected_umi = umi_map_for_cb.get(umi)
    if corrected_umi is None:
        stats['umi_not_in_map_for_cb'] += 1
        return None

    barcode_length = cb_length + umi_length

    # Build unaligned BAM records
    r1_bam = pysam.AlignedSegment()
    r1_bam.query_name = r1_name
    r1_bam.query_sequence = remaining_r1_seq
    # Biopython stores qualities as list[int]; slice past CB+UMI
    r1_bam.query_qualities = r1_record.letter_annotations["phred_quality"][barcode_length:]
    r1_bam.flag = 77                 # 0x004D: paired, first in pair, unmapped, mate unmapped
    r1_bam.reference_id = -1
    r1_bam.reference_start = -1
    r1_bam.mapping_quality = 0
    r1_bam.set_tag("CB", corrected_cb)
    r1_bam.set_tag("UB", corrected_umi)

    r2_bam = pysam.AlignedSegment()
    r2_bam.query_name = r2_name
    r2_bam.query_sequence = str(r2_record.seq)
    r2_bam.query_qualities = r2_record.letter_annotations["phred_quality"]
    r2_bam.flag = 141                # 0x008D: paired, second in pair, unmapped, mate unmapped
    r2_bam.reference_id = -1
    r2_bam.reference_start = -1
    r2_bam.mapping_quality = 0
    r2_bam.set_tag("CB", corrected_cb)
    r2_bam.set_tag("UB", corrected_umi)

    return r1_bam, r2_bam


def fastq_to_bam_with_barcodes(r1_fastq, r2_fastq, cb_umi_mapping_file, output_bam, num_cores=1, cb_length=16, umi_length=12):
    """
    Convert paired FASTQ files to unaligned BAM with CB/UB tags using multiple threads.

    Args:
        r1_fastq: path to R1 FASTQ(.gz)
        r2_fastq: path to R2 FASTQ(.gz)
        cb_umi_mapping_file: TSV(.gz) with 4 columns (rawCB, corrCB, rawUMI, corrUMI)
        output_bam: path to output BAM
        num_cores: threads
        cb_length, umi_length: lengths for parsing R1
    """
    print("Loading CB/UMI mapping...")
    cb_map, umi_map_by_raw_cb = load_cb_umi_mapping(cb_umi_mapping_file)

    stats = defaultdict(int)

    r1_open_func = gzip.open if r1_fastq.endswith('.gz') else open
    r2_open_func = gzip.open if r2_fastq.endswith('.gz') else open
    r1_mode = 'rt' if r1_fastq.endswith('.gz') else 'r'
    r2_mode = 'rt' if r2_fastq.endswith('.gz') else 'r'

    header = {
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'PG': [{'ID': 'nimble-fastq-to-bam', 'PN': 'nimble', 'VN': '1.1', 'CL': 'single-tsv, cb-scoped-umi'}]
    }

    print(f"Processing paired FASTQ files with {num_cores} threads...")

    try:
        with r1_open_func(r1_fastq, r1_mode) as r1_handle, \
             r2_open_func(r2_fastq, r2_mode) as r2_handle, \
             pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:

            r1_iter = SeqIO.parse(r1_handle, "fastq")
            r2_iter = SeqIO.parse(r2_handle, "fastq")

            with ThreadPoolExecutor(max_workers=num_cores) as executor:
                futures = {}

                for idx, (r1_record, r2_record) in enumerate(zip(r1_iter, r2_iter), start=1):
                    stats['total_pairs'] += 1
                    fut = executor.submit(process_pair, r1_record, r2_record,
                                          cb_map, umi_map_by_raw_cb, stats, cb_length, umi_length)
                    futures[fut] = True

                    # Throttle in-flight futures to limit memory use
                    if len(futures) >= num_cores * 100:
                        done_any = 0
                        for done in as_completed(list(futures)[:num_cores * 10]):
                            res = done.result()
                            if res:
                                r1_bam, r2_bam = res
                                bam_out.write(r1_bam)
                                bam_out.write(r2_bam)
                                stats['written_pairs'] += 1
                            del futures[done]
                            done_any += 1
                            if done_any >= num_cores * 10:
                                break

                    if stats['total_pairs'] % 1_000_000 == 0:
                        print(f"Processed {stats['total_pairs']} read pairs...")

                # Drain remaining futures
                for done in as_completed(list(futures.keys())):
                    res = done.result()
                    if res:
                        r1_bam, r2_bam = res
                        bam_out.write(r1_bam)
                        bam_out.write(r2_bam)
                        stats['written_pairs'] += 1
                    del futures[done]

    except Exception as e:
        print(f"Error during processing: {e}", file=sys.stderr)
        sys.exit(1)

    # Print summary
    print("\n=== Processing Statistics ===")
    ordered = [
        'total_pairs', 'written_pairs',
        'name_mismatch', 'too_short', 'no_remaining_seq',
        'cb_not_in_map', 'cb_has_no_umi_map', 'umi_not_in_map_for_cb'
    ]
    for k in ordered:
        print(f"{k.replace('_', ' ').capitalize()}: {stats.get(k, 0)}")
    # Print any other counters encountered
    for k, v in stats.items():
        if k not in ordered:
            print(f"{k.replace('_', ' ').capitalize()}: {v}")
    print(f"Output BAM written to: {output_bam}")