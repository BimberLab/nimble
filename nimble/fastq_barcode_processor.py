#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import islice, chain
import pysam
from Bio import SeqIO


def hamming_distance(s1, s2):
    """Calculate Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        return float("inf")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def build_hamming_index(whitelist):
    """
    Build an index mapping each valid CB to all possible 1-edit variants.
    Returns: dict[str variant] = set[str valid_cb]
    """
    bases = ["A", "C", "G", "T", "N"]
    hamming_index = defaultdict(set)

    for valid_cb in whitelist:
        for i in range(len(valid_cb)):
            for base in bases:
                if base != valid_cb[i]:
                    variant = valid_cb[:i] + base + valid_cb[i + 1 :]
                    hamming_index[variant].add(valid_cb)

    return hamming_index


def load_cb_whitelist(whitelist_path):
    """
    Load a cell barcode whitelist (one CB per line) and build a Hamming distance = 1 index.
    Returns: (whitelist_set, hamming_index)
    """
    whitelist = set()

    open_func = gzip.open if whitelist_path.endswith(".gz") else open
    mode = "rt" if whitelist_path.endswith(".gz") else "r"

    try:
        with open_func(whitelist_path, mode) as f:
            for line in f:
                line = line.strip()
                if line:
                    whitelist.add(line)
    except Exception as e:
        print(f"Error loading CB whitelist from {whitelist_path}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded whitelist from {whitelist_path}")
    print(f"  Valid cell barcodes: {len(whitelist)}")

    print("Building Hamming distance index...")
    hamming_index = build_hamming_index(whitelist)
    print(f"  Indexed {len(hamming_index)} variants")

    return whitelist, hamming_index


def correct_cell_barcode(raw_cb, quality_scores, whitelist, hamming_index, correction_cache):
    """
    Correct a raw cell barcode using 10x-style correction:
      1) perfect match
      2) candidates at Hamming distance = 1 (via index)
      3) if multiple, choose the candidate whose differing position has the lowest Q
    Returns corrected CB or None.
    """
    if raw_cb in correction_cache:
        return correction_cache[raw_cb]

    if raw_cb in whitelist:
        correction_cache[raw_cb] = raw_cb
        return raw_cb

    candidates = hamming_index.get(raw_cb, set())
    if not candidates:
        correction_cache[raw_cb] = None
        return None

    if len(candidates) == 1:
        corrected = next(iter(candidates))
        correction_cache[raw_cb] = corrected
        return corrected

    best_candidate = None
    lowest_quality = float("inf")

    for candidate in candidates:
        for i, (raw_base, cand_base) in enumerate(zip(raw_cb, candidate)):
            if raw_base != cand_base:
                qual = quality_scores[i]
                if qual < lowest_quality:
                    lowest_quality = qual
                    best_candidate = candidate
                break

    correction_cache[raw_cb] = best_candidate
    return best_candidate

def infer_umi_and_tso_lengths_from_r1_records(
    r1_records,
    search_string="TTTCTTATATGGG",
    cb_length=16,
    min_records=10,
):
    """
    Infer UMI length by locating the TSO motif in R1 reads.
    Assumes structure: [CB][UMI][TSO...][cDNA...]

    Resolution strategy:
      - Ignore reads with no TSO hit (pos == -1).
      - Allow multiple observed TSO positions; pick the *mode* (majority) position.
      - If fewer than min_records hits are available (e.g. EOF), use as many as we have.
    """

    tso_starts = []
    for rec in r1_records:
        seq = str(rec.seq)
        pos = seq.find(search_string)
        if pos == -1:
            continue  # skip non-hits for inference
        tso_starts.append(pos)

    if len(tso_starts) == 0:
        raise ValueError(
            f"Unable to infer UMI length: no reads contained TSO motif {search_string}"
        )

    counts = Counter(tso_starts)
    # Choose the most common; break ties by choosing the smallest position
    mode_pos, mode_count = min(
        (pos, cnt) for pos, cnt in counts.items() if cnt == max(counts.values())
    )
    total_hits = len(tso_starts)
    top = counts.most_common(10)
    top_str = ", ".join([f"{p}:{c}" for p, c in top])
    print(f"TSO hit count used for inference: {total_hits}")
    print(f"TSO positions (top): {top_str}")
    if total_hits < min_records:
        print(
            f"WARNING: only {total_hits} reads contained the TSO motif; "
            f"min_records={min_records}. Proceeding with majority position."
        )

    umi_length = mode_pos - cb_length
    if umi_length <= 0:
        raise ValueError(
            f"Inferred UMI length <= 0 (mode_tso_start={mode_pos}, cb_length={cb_length}). "
            f"Check motif/structure."
        )

    tso_length = len(search_string)
    print(f"Inferred UMI length: {umi_length} (from majority TSO start {mode_pos})")
    print(f"Using TSO length: {tso_length} (len of search_string)")
    return umi_length, tso_length

def parse_10x_barcode_from_r1(sequence, cb_length=16, umi_length=12, tso_length=0):
    """
    Parse CB and UMI from R1 sequence, then drop optional TSO.
    Returns (cell_barcode, umi, remaining_sequence).
    """
    prefix_len = cb_length + umi_length + tso_length
    if len(sequence) < prefix_len:
        return None, None, ""
    cell_barcode = sequence[:cb_length]
    umi = sequence[cb_length : cb_length + umi_length]
    remaining_sequence = sequence[prefix_len:]
    return cell_barcode, umi, remaining_sequence


def process_pair(
    r1_record,
    r2_record,
    whitelist,
    hamming_index,
    correction_cache,
    stats,
    cb_length=16,
    umi_length=12,
    tso_length=0,
):
    """
    Process a single FASTQ pair; return (r1_bam, r2_bam) or None if skipped.
    - Correct CB (Hamming distance 0/1) using CB qualities.
    - Use raw UMI.
    - Drop (CB+UMI+TSO) from R1 sequence and qualities consistently.
    """
    r1_name = r1_record.id.removesuffix("/1")
    r2_name = r2_record.id.removesuffix("/2")
    if r1_name != r2_name:
        stats["name_mismatch"] += 1
        return None

    r1_seq = str(r1_record.seq)
    raw_cb, umi, remaining_r1_seq = parse_10x_barcode_from_r1(
        r1_seq, cb_length=cb_length, umi_length=umi_length, tso_length=tso_length
    )
    if raw_cb is None or umi is None:
        stats["too_short"] += 1
        return None
    if len(remaining_r1_seq) == 0:
        stats["no_remaining_seq"] += 1
        return None

    # CB qualities: first cb_length bases
    cb_quality_scores = r1_record.letter_annotations["phred_quality"][:cb_length]

    corrected_cb = correct_cell_barcode(raw_cb, cb_quality_scores, whitelist, hamming_index, correction_cache)
    if corrected_cb is None:
        stats["cb_no_correction"] += 1
        return None

    if corrected_cb == raw_cb:
        stats["cb_perfect_match"] += 1
    else:
        stats["cb_corrected"] += 1

    prefix_len = cb_length + umi_length + tso_length

    r1_bam = pysam.AlignedSegment()
    r1_bam.query_name = r1_name
    r1_bam.query_sequence = remaining_r1_seq
    r1_bam.query_qualities = r1_record.letter_annotations["phred_quality"][prefix_len:]
    r1_bam.flag = 77  # paired, first in pair, unmapped, mate unmapped
    r1_bam.reference_id = -1
    r1_bam.reference_start = -1
    r1_bam.mapping_quality = 0
    r1_bam.set_tag("CB", corrected_cb)
    r1_bam.set_tag("UB", umi)

    r2_bam = pysam.AlignedSegment()
    r2_bam.query_name = r2_name
    r2_bam.query_sequence = str(r2_record.seq)
    r2_bam.query_qualities = r2_record.letter_annotations["phred_quality"]
    r2_bam.flag = 141  # paired, second in pair, unmapped, mate unmapped
    r2_bam.reference_id = -1
    r2_bam.reference_start = -1
    r2_bam.mapping_quality = 0
    r2_bam.set_tag("CB", corrected_cb)
    r2_bam.set_tag("UB", umi)

    return r1_bam, r2_bam


def fastq_to_bam_with_barcodes(
    r1_fastq,
    r2_fastq,
    cb_whitelist_file,
    output_bam,
    num_cores=1,
    cb_length=16,
    umi_length=None,
    infer_umi=True,
    tso_search_string="TTTCTTATATGGG",
    infer_prefix_pairs=200,
    min_records_with_tso=10,
):
    """
    Convert paired FASTQ files to unaligned BAM with CB/UB tags using multiple threads.
    - Corrects CB using whitelist-based Hamming distance 0/1 correction with CB qualities.
    - Optionally infers UMI length by locating a TSO motif in R1 and removes the TSO from R1.
    
    Args:
        r1_fastq, r2_fastq: FASTQ(.gz) paths
        cb_whitelist_file: whitelist file
        output_bam: output BAM path
        num_cores: threads
        cb_length: cell barcode length
        cb_whitelist_file: whitelist file
        output_bam: output BAM path
        num_cores: threads
        cb_length: cell barcode length
        umi_length: if infer_umi=False, required
        infer_umi: infer umi_length and remove tso motif
        tso_search_string: TSO motif used for inference; tso_length = len(tso_search_string)
        infer_prefix_pairs: number of pairs to buffer for inference
        min_records_with_tso: minimum reads containing motif to accept inference
    """
    print("Loading cell barcode whitelist...")
    whitelist, hamming_index = load_cb_whitelist(cb_whitelist_file)

    correction_cache = {}
    stats = defaultdict(int)

    r1_open_func = gzip.open if r1_fastq.endswith(".gz") else open
    r2_open_func = gzip.open if r2_fastq.endswith(".gz") else open
    r1_mode = "rt" if r1_fastq.endswith(".gz") else "r"
    r2_mode = "rt" if r2_fastq.endswith(".gz") else "r"

    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "PG": [
            {
                "ID": "nimble-fastq-to-bam",
                "PN": "nimble",
                "VN": "1.2",
                "CL": "whitelist-based CB correction (+ optional UMI inference + TSO removal)",
            }
        ],
    }

    print(f"Processing paired FASTQ files with {num_cores} threads...")

    try:
        with r1_open_func(r1_fastq, r1_mode) as r1_handle, \
             r2_open_func(r2_fastq, r2_mode) as r2_handle, \
             pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:

            r1_iter = SeqIO.parse(r1_handle, "fastq")
            r2_iter = SeqIO.parse(r2_handle, "fastq")

            prefix_pairs = []
            tso_hit_r1 = []

            target_hits = infer_prefix_pairs if infer_umi else 0

            while True:
                try:
                    r1_record = next(r1_iter)
                    r2_record = next(r2_iter)
                except StopIteration:
                    break

                prefix_pairs.append((r1_record, r2_record))

                if infer_umi:
                    pos = str(r1_record.seq).find(tso_search_string)
                    if pos != -1:
                        tso_hit_r1.append(r1_record)

                    # Stop once we have enough usable reads for inference
                    if len(tso_hit_r1) >= target_hits:
                        break

                # If not inferring, we only need a minimal buffer to avoid empty input
                if not infer_umi and len(prefix_pairs) >= 1:
                    break

            if not prefix_pairs:
                raise ValueError("Input FASTQs are empty or could not be read")

            if infer_umi:
                inferred_umi_length, tso_length = infer_umi_and_tso_lengths_from_r1_records(
                    tso_hit_r1,
                    search_string=tso_search_string,
                    cb_length=cb_length,
                    min_records=min_records_with_tso,
                )
                umi_length_use = inferred_umi_length
            else:
                if umi_length is None:
                    raise ValueError("umi_length must be provided when infer_umi=False")
                umi_length_use = int(umi_length)
                tso_length = 0

            r1_full = chain((p[0] for p in prefix_pairs), r1_iter)
            r2_full = chain((p[1] for p in prefix_pairs), r2_iter)

            with ThreadPoolExecutor(max_workers=num_cores) as executor:
                futures = {}

                for idx, (r1_record, r2_record) in enumerate(zip(r1_full, r2_full), start=1):
                    stats["total_pairs"] += 1
                    fut = executor.submit(
                        process_pair,
                        r1_record,
                        r2_record,
                        whitelist,
                        hamming_index,
                        correction_cache,
                        stats,
                        cb_length,
                        umi_length_use,
                        tso_length,
                    )
                    futures[fut] = True

                    if len(futures) >= num_cores * 100:
                        done_any = 0
                        for done in as_completed(list(futures)[: num_cores * 10]):
                            res = done.result()
                            if res:
                                r1_bam, r2_bam = res
                                bam_out.write(r1_bam)
                                bam_out.write(r2_bam)
                                stats["written_pairs"] += 1
                            del futures[done]
                            done_any += 1
                            if done_any >= num_cores * 10:
                                break

                    if stats["total_pairs"] % 1_000_000 == 0:
                        print(f"Processed {stats['total_pairs']} read pairs...")

                for done in as_completed(list(futures.keys())):
                    res = done.result()
                    if res:
                        r1_bam, r2_bam = res
                        bam_out.write(r1_bam)
                        bam_out.write(r2_bam)
                        stats["written_pairs"] += 1
                    del futures[done]

    except Exception as e:
        print(f"Error during processing: {e}", file=sys.stderr)
        sys.exit(1)

    print("\n=== Processing Statistics ===")
    print(f"Total read pairs: {stats.get('total_pairs', 0)}")
    print(f"Written pairs: {stats.get('written_pairs', 0)}")

    print("\nCell Barcode Correction:")
    print(f"  Perfect matches: {stats.get('cb_perfect_match', 0)}")
    print(f"  Corrected (1-edit): {stats.get('cb_corrected', 0)}")
    print(f"  No valid correction: {stats.get('cb_no_correction', 0)}")

    total_cb_processed = (
        stats.get("cb_perfect_match", 0)
        + stats.get("cb_corrected", 0)
        + stats.get("cb_no_correction", 0)
    )

    if total_cb_processed > 0:
        perfect_pct = 100.0 * stats.get("cb_perfect_match", 0) / total_cb_processed
        corrected_pct = 100.0 * stats.get("cb_corrected", 0) / total_cb_processed
        dropped_pct = 100.0 * stats.get("cb_no_correction", 0) / total_cb_processed
        print(f"  Correction rate: {perfect_pct:.2f}% perfect, {corrected_pct:.2f}% corrected, {dropped_pct:.2f}% dropped")

        if dropped_pct == 100.0:
            raise ValueError("There were no passing cell barcodes. This likely indicates an error with the whitelist.")
    else:
        raise ValueError("No cell barcodes were processed. There is likely a problem with the input files.")

    print("\nOther filters:")
    print(f"  Name mismatch: {stats.get('name_mismatch', 0)}")
    print(f"  Too short: {stats.get('too_short', 0)}")
    print(f"  No remaining sequence: {stats.get('no_remaining_seq', 0)}")

    print(f"\nCorrection cache size: {len(correction_cache)} unique raw CBs")
    print(f"\nOutput BAM written to: {output_bam}")