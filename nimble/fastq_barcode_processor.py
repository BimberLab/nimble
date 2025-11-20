#!/usr/bin/env python3
import gzip
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import pysam
from Bio import SeqIO


def hamming_distance(s1, s2):
    """Calculate Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def build_hamming_index(whitelist):
    """
    Build an index mapping each valid CB to all possible 1-edit variants.
    This allows efficient lookup of correction candidates.
    
    Returns:
        dict[str variant] = set[str valid_cb]
    """
    bases = ['A', 'C', 'G', 'T', 'N']
    hamming_index = defaultdict(set)
    
    for valid_cb in whitelist:
        # For each position, generate all possible single-base substitutions
        for i in range(len(valid_cb)):
            for base in bases:
                if base != valid_cb[i]:
                    variant = valid_cb[:i] + base + valid_cb[i+1:]
                    hamming_index[variant].add(valid_cb)
    
    return hamming_index


def load_cb_whitelist(whitelist_path):
    """
    Load a cell barcode whitelist (one CB per line).
    Build a Hamming distance = 1 index for efficient correction.

    Returns:
        whitelist: set of valid cell barcodes
        hamming_index: dict mapping variants to valid CBs
    """
    whitelist = set()
    
    open_func = gzip.open if whitelist_path.endswith('.gz') else open
    mode = 'rt' if whitelist_path.endswith('.gz') else 'r'

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
    1. Check for perfect match
    2. If not, find all candidates with Hamming distance = 1
    3. Among candidates, select the one where the differing base has the lowest quality score
    
    Args:
        raw_cb: Raw cell barcode string
        quality_scores: List of quality scores for the CB region
        whitelist: Set of valid cell barcodes
        hamming_index: Dict mapping variants to valid CBs
        correction_cache: Dict for caching corrections
    
    Returns:
        Corrected CB string, or None if no valid correction found
    """
    # Check cache first
    if raw_cb in correction_cache:
        return correction_cache[raw_cb]
    
    # Check for perfect match
    if raw_cb in whitelist:
        correction_cache[raw_cb] = raw_cb
        return raw_cb
    
    # Find candidates with Hamming distance = 1
    candidates = hamming_index.get(raw_cb, set())
    
    if not candidates:
        correction_cache[raw_cb] = None
        return None
    
    if len(candidates) == 1:
        # Only one candidate, use it
        corrected = next(iter(candidates))
        correction_cache[raw_cb] = corrected
        return corrected
    
    # Multiple candidates: select based on quality scores
    # Find the differing position and pick candidate with lowest quality at that position
    best_candidate = None
    lowest_quality = float('inf')
    
    for candidate in candidates:
        # Find the position where they differ
        for i, (raw_base, cand_base) in enumerate(zip(raw_cb, candidate)):
            if raw_base != cand_base:
                qual = quality_scores[i]
                if qual < lowest_quality:
                    lowest_quality = qual
                    best_candidate = candidate
                break
    
    correction_cache[raw_cb] = best_candidate
    return best_candidate


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


def process_pair(r1_record, r2_record, whitelist, hamming_index, correction_cache, stats, cb_length=16, umi_length=12):
    """
    Process a single FASTQ pair and return (r1_bam, r2_bam), or None if skipped.
    Correction logic:
      - Correct CB using 10x-style Hamming distance = 1 correction with quality scores
      - Use raw UMI (no correction)
    """
    # Normalize names (handle /1 and /2 suffixes if present)
    r1_name = r1_record.id.removesuffix('/1')
    r2_name = r2_record.id.removesuffix('/2')
    if r1_name != r2_name:
        stats['name_mismatch'] += 1
        return None

    r1_seq = str(r1_record.seq)
    raw_cb, umi, remaining_r1_seq = parse_10x_barcode_from_r1(r1_seq, cb_length, umi_length)
    if raw_cb is None or umi is None:
        stats['too_short'] += 1
        return None
    if len(remaining_r1_seq) == 0:
        stats['no_remaining_seq'] += 1
        return None

    # Get quality scores for the CB region
    cb_quality_scores = r1_record.letter_annotations["phred_quality"][:cb_length]
    
    # Correct cell barcode using 10x-style correction
    corrected_cb = correct_cell_barcode(raw_cb, cb_quality_scores, whitelist, hamming_index, correction_cache)
    
    if corrected_cb is None:
        stats['cb_no_correction'] += 1
        return None
    
    # Track correction statistics
    if corrected_cb == raw_cb:
        stats['cb_perfect_match'] += 1
    else:
        stats['cb_corrected'] += 1

    barcode_length = cb_length + umi_length

    # Build unaligned BAM records with corrected CB and raw UMI
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
    r1_bam.set_tag("UB", umi)  # Use raw UMI

    r2_bam = pysam.AlignedSegment()
    r2_bam.query_name = r2_name
    r2_bam.query_sequence = str(r2_record.seq)
    r2_bam.query_qualities = r2_record.letter_annotations["phred_quality"]
    r2_bam.flag = 141                # 0x008D: paired, second in pair, unmapped, mate unmapped
    r2_bam.reference_id = -1
    r2_bam.reference_start = -1
    r2_bam.mapping_quality = 0
    r2_bam.set_tag("CB", corrected_cb)
    r2_bam.set_tag("UB", umi)  # Use raw UMI

    return r1_bam, r2_bam


def fastq_to_bam_with_barcodes(r1_fastq, r2_fastq, cb_whitelist_file, output_bam, num_cores=1, cb_length=16, umi_length=12):
    """
    Convert paired FASTQ files to unaligned BAM with CB/UB tags using multiple threads.
    Performs 10x-style cell barcode correction using a whitelist.

    Args:
        r1_fastq: path to R1 FASTQ(.gz)
        r2_fastq: path to R2 FASTQ(.gz)
        cb_whitelist_file: path to cell barcode whitelist file (one CB per line, .gz or plain text)
        output_bam: path to output BAM
        num_cores: threads
        cb_length: length of cell barcode (default 16)
        umi_length: length of UMI (default 12)
    """
    print("Loading cell barcode whitelist...")
    whitelist, hamming_index = load_cb_whitelist(cb_whitelist_file)
    
    # Initialize correction cache (shared across threads, but thread-safe via GIL for dict operations)
    correction_cache = {}
    
    stats = defaultdict(int)

    r1_open_func = gzip.open if r1_fastq.endswith('.gz') else open
    r2_open_func = gzip.open if r2_fastq.endswith('.gz') else open
    r1_mode = 'rt' if r1_fastq.endswith('.gz') else 'r'
    r2_mode = 'rt' if r2_fastq.endswith('.gz') else 'r'

    header = {
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'PG': [{'ID': 'nimble-fastq-to-bam', 'PN': 'nimble', 'VN': '1.2', 'CL': 'whitelist-based CB correction'}]
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
                                          whitelist, hamming_index, correction_cache, stats, 
                                          cb_length, umi_length)
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

    # Print summary with correction statistics
    print("\n=== Processing Statistics ===")
    print(f"Total read pairs: {stats.get('total_pairs', 0)}")
    print(f"Written pairs: {stats.get('written_pairs', 0)}")
    print(f"\nCell Barcode Correction:")
    print(f"  Perfect matches: {stats.get('cb_perfect_match', 0)}")
    print(f"  Corrected (1-edit): {stats.get('cb_corrected', 0)}")
    print(f"  No valid correction: {stats.get('cb_no_correction', 0)}")
    
    total_cb_processed = stats.get('cb_perfect_match', 0) + stats.get('cb_corrected', 0) + stats.get('cb_no_correction', 0)
    if total_cb_processed > 0:
        perfect_pct = 100.0 * stats.get('cb_perfect_match', 0) / total_cb_processed
        corrected_pct = 100.0 * stats.get('cb_corrected', 0) / total_cb_processed
        dropped_pct = 100.0 * stats.get('cb_no_correction', 0) / total_cb_processed
        print(f"  Correction rate: {perfect_pct:.2f}% perfect, {corrected_pct:.2f}% corrected, {dropped_pct:.2f}% dropped")
    
    print(f"\nOther filters:")
    print(f"  Name mismatch: {stats.get('name_mismatch', 0)}")
    print(f"  Too short: {stats.get('too_short', 0)}")
    print(f"  No remaining sequence: {stats.get('no_remaining_seq', 0)}")
    
    # Print unique cache size
    print(f"\nCorrection cache size: {len(correction_cache)} unique raw CBs")
    
    print(f"\nOutput BAM written to: {output_bam}")
