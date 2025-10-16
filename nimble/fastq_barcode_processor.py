#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import pysam
from Bio import SeqIO

def load_barcode_mapping(tsv_path):
    """
    Load a TSV file with Raw->Corrected barcode mappings into a dictionary.
    Handles gzipped files automatically.
    """
    barcode_map = {}

    open_func = gzip.open if tsv_path.endswith('.gz') else open
    mode = 'rt' if tsv_path.endswith('.gz') else 'r'

    try:
        with open_func(tsv_path, mode) as f:
            header = f.readline().strip().split('\t')
            if len(header) != 2:
                raise ValueError(f"Expected 2 columns in TSV header, got {len(header)}")

            for line_num, line in enumerate(f, start=2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    print(f"Warning: Skipping malformed line {line_num} in {tsv_path}: {line}")
                    continue
                raw_barcode, corrected_barcode = parts
                barcode_map[raw_barcode] = corrected_barcode

    except Exception as e:
        print(f"Error loading barcode mapping from {tsv_path}: {e}")
        sys.exit(1)

    print(f"Loaded {len(barcode_map)} barcode mappings from {tsv_path}")
    return barcode_map


def parse_10x_barcode_from_r1(sequence, cb_length=16, umi_length=12):
    """
    Parse cell barcode and UMI from 10x R1 read sequence.
    """
    if len(sequence) < cb_length + umi_length:
        return None, None, ""

    cell_barcode = sequence[:cb_length]
    umi = sequence[cb_length:cb_length + umi_length]
    remaining_sequence = sequence[cb_length + umi_length:]

    return cell_barcode, umi, remaining_sequence


def process_pair(r1_record, r2_record, cb_map, umi_map):
    """
    Process a single FASTQ pair and return (r1_bam, r2_bam), or None if skipped.
    """
    r1_name = r1_record.id.removesuffix('/1')
    r2_name = r2_record.id.removesuffix('/2')
    if r1_name != r2_name:
        return None

    r1_seq = str(r1_record.seq)
    cell_barcode, umi, remaining_r1_seq = parse_10x_barcode_from_r1(r1_seq)
    if cell_barcode is None or umi is None or len(remaining_r1_seq) == 0:
        return None

    corrected_cb = cb_map.get(cell_barcode)
    corrected_umi = umi_map.get(umi)
    if corrected_cb is None or corrected_umi is None:
        return None

    barcode_length = len(cell_barcode) + len(umi)

    # R1 BAM
    r1_bam = pysam.AlignedSegment()
    r1_bam.query_name = r1_name
    r1_bam.query_sequence = remaining_r1_seq
    r1_bam.query_qualities = r1_record.letter_annotations["phred_quality"][barcode_length:]
    r1_bam.flag = 77
    r1_bam.reference_id = -1
    r1_bam.reference_start = -1
    r1_bam.mapping_quality = 0
    r1_bam.set_tag("CB", corrected_cb)
    r1_bam.set_tag("UB", corrected_umi)

    # R2 BAM
    r2_bam = pysam.AlignedSegment()
    r2_bam.query_name = r2_name
    r2_bam.query_sequence = str(r2_record.seq)
    r2_bam.query_qualities = r2_record.letter_annotations["phred_quality"]
    r2_bam.flag = 141
    r2_bam.reference_id = -1
    r2_bam.reference_start = -1
    r2_bam.mapping_quality = 0
    r2_bam.set_tag("CB", corrected_cb)
    r2_bam.set_tag("UB", corrected_umi)

    return r1_bam, r2_bam


def fastq_to_bam_with_barcodes(r1_fastq, r2_fastq, cb_mapping_file, umi_mapping_file, output_bam, num_cores=1):
    """
    Convert paired FASTQ files to unaligned BAM with CB/UB tags using multiple threads.
    """
    print("Loading barcode mappings...")
    cb_map = load_barcode_mapping(cb_mapping_file)
    umi_map = load_barcode_mapping(umi_mapping_file)

    stats = defaultdict(int)

    r1_open_func = gzip.open if r1_fastq.endswith('.gz') else open
    r2_open_func = gzip.open if r2_fastq.endswith('.gz') else open
    r1_mode = 'rt' if r1_fastq.endswith('.gz') else 'r'
    r2_mode = 'rt' if r2_fastq.endswith('.gz') else 'r'

    header = {
        'HD': {'VN': '1.4', 'SO': 'queryname'},
        'PG': [{'ID': 'nimble-fastq-to-bam', 'PN': 'nimble', 'VN': '1.0'}]
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
                    stats['total_reads'] += 1
                    fut = executor.submit(process_pair, r1_record, r2_record, cb_map, umi_map)
                    futures[fut] = idx

                    # Throttle in-flight futures to limit memory use
                    if len(futures) >= num_cores * 100:
                        for done in as_completed(list(futures)[:num_cores]):
                            result = done.result()
                            if result:
                                r1_bam, r2_bam = result
                                bam_out.write(r1_bam)
                                bam_out.write(r2_bam)
                                stats['written_pairs'] += 1
                            del futures[done]

                    if stats['total_reads'] % 1000000 == 0:
                        print(f"Processed {stats['total_reads']} reads...")

                # Drain remaining futures
                for done in as_completed(futures):
                    result = done.result()
                    if result:
                        r1_bam, r2_bam = result
                        bam_out.write(r1_bam)
                        bam_out.write(r2_bam)
                        stats['written_pairs'] += 1

    except Exception as e:
        print(f"Error during processing: {e}")
        sys.exit(1)

    # Print summary
    print("\n=== Processing Statistics ===")
    for k, v in stats.items():
        print(f"{k.replace('_', ' ').capitalize()}: {v}")
    print(f"Output BAM written to: {output_bam}")