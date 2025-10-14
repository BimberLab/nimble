#!/usr/bin/env python
import gzip
import os
import sys
from collections import defaultdict
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def load_barcode_mapping(tsv_path):
    """
    Load a TSV file with Raw->Corrected barcode mappings into a dictionary.
    Automatically handles gzipped files based on .gz extension.
    
    Args:
        tsv_path: Path to TSV file (can be .txt or .txt.gz)
        
    Returns:
        dict: Dictionary mapping raw barcodes to corrected barcodes
    """
    barcode_map = {}
    
    # Determine if file is gzipped based on extension
    if tsv_path.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
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
    Standard 10x format: [16bp CB][10-12bp UMI][remaining sequence]
    
    Args:
        sequence: R1 sequence string
        cb_length: Length of cell barcode (default 16bp for 10x)
        umi_length: Length of UMI (default 12bp, but can be 10)
        
    Returns:
        tuple: (cell_barcode, umi, remaining_sequence)
    """
    if len(sequence) < cb_length + umi_length:
        return None, None, ""
    
    cell_barcode = sequence[:cb_length]
    umi = sequence[cb_length:cb_length + umi_length]
    remaining_sequence = sequence[cb_length + umi_length:]
    
    return cell_barcode, umi, remaining_sequence


def fastq_to_bam_with_barcodes(r1_fastq, r2_fastq, cb_mapping_file, umi_mapping_file, output_bam, num_cores=1):
    """
    Convert paired FASTQ files with separate barcode mappings to unaligned BAM with CB/UB tags.
    
    Args:
        r1_fastq: Path to R1 FASTQ file
        r2_fastq: Path to R2 FASTQ file  
        cb_mapping_file: Path to cell barcode mapping TSV
        umi_mapping_file: Path to UMI mapping TSV
        output_bam: Path for output BAM file
        num_cores: Number of cores for processing
    """
    
    print("Loading barcode mappings...")
    cb_map = load_barcode_mapping(cb_mapping_file)
    umi_map = load_barcode_mapping(umi_mapping_file)
    
    print(f"Processing paired FASTQ files: {r1_fastq}, {r2_fastq}")
    
    # Statistics tracking
    stats = {
        'total_reads': 0,
        'skipped_short_r1': 0,
        'missing_barcode': 0,
        'missing_cb_mapping': 0,
        'missing_umi_mapping': 0,
        'written_pairs': 0
    }
    
    # Determine if FASTQ files are gzipped
    r1_open_func = gzip.open if r1_fastq.endswith('.gz') else open
    r2_open_func = gzip.open if r2_fastq.endswith('.gz') else open
    r1_mode = 'rt' if r1_fastq.endswith('.gz') else 'r'
    r2_mode = 'rt' if r2_fastq.endswith('.gz') else 'r'
    
    # Create BAM header - minimal header for unaligned reads
    header = {
        'HD': {'VN': '1.4', 'SO': 'queryname'},
        'PG': [{'ID': 'nimble-fastq-to-bam', 'PN': 'nimble', 'VN': '1.0'}]
    }
    
    try:
        with r1_open_func(r1_fastq, r1_mode) as r1_handle, \
             r2_open_func(r2_fastq, r2_mode) as r2_handle, \
             pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
            
            # Use SeqIO to parse FASTQ files simultaneously
            r1_records = SeqIO.parse(r1_handle, "fastq")
            r2_records = SeqIO.parse(r2_handle, "fastq") 
            
            for r1_record, r2_record in zip(r1_records, r2_records):
                stats['total_reads'] += 1
                
                # Ensure read names match (without /1 /2 suffixes)
                r1_name = r1_record.id.removesuffix('/1')
                r2_name = r2_record.id.removesuffix('/2')
                if r1_name != r2_name:
                    print(f"Warning: Skipping pair, read name mismatch: {r1_name} vs {r2_name}")
                    continue
                
                # Parse barcode information from R1
                r1_seq = str(r1_record.seq)
                cell_barcode, umi, remaining_r1_seq = parse_10x_barcode_from_r1(r1_seq)
                
                if cell_barcode is None or umi is None:
                    stats['missing_barcode'] += 1
                    print(f"Warning: Skipping pair, missing barcode(s) for R1: {r1_name}")
                    continue
                
                # Skip if no sequence remains in R1 after barcode stripping
                if len(remaining_r1_seq) == 0:
                    stats['skipped_short_r1'] += 1
                    print(f"Warning: Skipping pair, short R1 read after barcode trimming: {r1_name}")
                    continue
                
                # Calculate barcode length for quality trimming
                barcode_length = len(cell_barcode) + len(umi)
                
                # Apply barcode corrections
                corrected_cb = cb_map.get(cell_barcode)
                corrected_umi = umi_map.get(umi)
                
                if corrected_cb is None:
                    stats['missing_cb_mapping'] += 1
                    print(f"Warning: Skipping pair, no mapping found for cell barcode: {cell_barcode} in R1: {r1_name}")
                    continue
                    
                if corrected_umi is None:
                    stats['missing_umi_mapping'] += 1
                    print(f"Warning: Skipping pair, no mapping found for UMI: {umi} in R1: {r1_name}")
                    continue

                # Create BAM records for R1 and R2
                r1_bam = pysam.AlignedSegment()
                r1_bam.query_name = r1_name
                r1_bam.query_sequence = remaining_r1_seq
                r1_bam.query_qualities = r1_record.letter_annotations["phred_quality"][barcode_length:]  # Skip barcode qualities
                r1_bam.flag = 77  # Read paired, read unmapped, mate unmapped, first in pair
                r1_bam.reference_id = -1
                r1_bam.reference_start = -1
                r1_bam.mapping_quality = 0
                r1_bam.set_tag("CB", corrected_cb)
                r1_bam.set_tag("UB", corrected_umi)
                
                r2_bam = pysam.AlignedSegment()
                r2_bam.query_name = r2_name  
                r2_bam.query_sequence = str(r2_record.seq)
                r2_bam.query_qualities = r2_record.letter_annotations["phred_quality"]
                r2_bam.flag = 141  # Read paired, read unmapped, mate unmapped, second in pair
                r2_bam.reference_id = -1
                r2_bam.reference_start = -1
                r2_bam.mapping_quality = 0
                r2_bam.set_tag("CB", corrected_cb)
                r2_bam.set_tag("UB", corrected_umi)
                
                # Write both records to BAM
                bam_out.write(r1_bam)
                bam_out.write(r2_bam)
                
                stats['written_pairs'] += 1
                
                # Progress reporting
                if stats['total_reads'] % 100000 == 0:
                    print(f"Processed {stats['total_reads']} reads, written {stats['written_pairs']} pairs")
    
    except Exception as e:
        print(f"Error processing FASTQ files: {e}")
        sys.exit(1)
    
    # Print final statistics
    print("\n=== Processing Statistics ===")
    print(f"Total reads processed: {stats['total_reads']}")
    print(f"Written read pairs: {stats['written_pairs']}")
    print(f"Skipped (short R1): {stats['skipped_short_r1']}")
    print(f"Missing CB mappings: {stats['missing_cb_mapping']}")
    print(f"Missing UMI mappings: {stats['missing_umi_mapping']}")
    
    if stats['missing_cb_mapping'] > 0 or stats['missing_umi_mapping'] > 0:
        print("Warning: Some reads were skipped due to missing barcode mappings")
    
    print(f"Output BAM written to: {output_bam}")
