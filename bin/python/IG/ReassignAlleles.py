#!/usr/bin/env python3

"""
Novel Allele Reassignment:
  - Reads AIRR‑formatted rearrangement records using read_rearrangement.
  - Reads novel allele info (TSV) and germline FASTA.
  - Subsets records whose call column (default "v_call") contains any novel gene.
  - Builds a germline_novel dictionary from the union of novel alleles and assignments.
  - For each selected record, computes the best candidate allele using Hamming distance.
    If the candidate is not already in the call (a comma‑separated string), it is merged
    into the call at the appropriate order (with lower distance first).
  - If a likelihood column is provided, inserts a None (as NA) at the candidate's insertion index.
  - Writes the updated records using derive_rearrangement.
"""

import os
import argparse
import re
import pandas as pd
from Bio import SeqIO
from airr import read_rearrangement, derive_rearrangement
import sys

# -- Constants and Defaults --
LIKELIHOOD_RE = re.compile(r'\s+', re.I)


# -- Utility Functions --

def hamming_distance(seq1, seq2):
    """Compute Hamming distance between two sequences (assumes equal length)."""
    return sum(1 for a, b in zip(seq1, seq2) if a != b)


def merge_assignment_with_index(orig_call, candidate, sample_seq, allele_db, max_length=3):
    """
    Merge the candidate allele into the original assignment string.
    
    The candidate is inserted into the existing comma‑separated call so that alleles with lower 
    Hamming distances (to the sample sequence) appear first.
    
    Returns a tuple (merged_call, insertion_index). If the candidate is already present, 
    returns (orig_call, None).
    """
    orig_list = [x.strip() for x in orig_call.split(",") if x.strip()]
    
    # If candidate contains multiple alleles (in a tie), take the first for ordering.
    candidate_alleles = [x.strip() for x in candidate.split(",") if x.strip()]
    if not candidate_alleles:
        return orig_call, None
    candidate_allele = candidate_alleles[0]
    
    if candidate_allele in orig_list:
        return orig_call, None

    candidate_distance = hamming_distance(sample_seq, allele_db.get(candidate_allele, sample_seq))
    
    # Compute distances for alleles already in the original call.
    distances = []
    for allele in orig_list:
        d = hamming_distance(sample_seq, allele_db.get(allele, sample_seq))
        distances.append(d)
    
    # Determine insertion index: count alleles with distance lower than candidate's.
    insertion_index = 0
    for d in distances:
        if d < candidate_distance:
            insertion_index += 1
        else:
            break
    
    new_list = orig_list.copy()
    new_list.insert(insertion_index, candidate_allele)
    if len(new_list) > max_length:
        new_list = new_list[:max_length]
    merged_call = ",".join(new_list)
    return merged_call, insertion_index


def reassign_novel_allele(records, germline, call="v_call", seq="sequence_alignment",
                          max_call=3, likelihood_col=None):
    """
    Reassigns the allele for each record (in-place) using the provided germline sequences.
    
    For each record in the list of record dicts:
      - Compute Hamming distances between the sample's sequence (record[seq]) and every allele in germline.
      - Select the best candidate allele (or join tied candidates).
      - If the candidate is not already in the record's call (a comma‑separated string), merge it into the call.
      - If a likelihood column is provided, insert a None at the same index.
    
    Returns the list of updated records.
    """
    for record in records:
        sample_seq = record.get(seq, "")
        orig_call = record.get(call, "")
        
        # Compute distances for every allele in germline.
        distances = {allele: hamming_distance(sample_seq, allele_seq) 
                    for allele, allele_seq in germline.items()}
        
        if not distances:
            continue
            
        min_distance = min(distances.values())
        best_alleles = [allele for allele, d in distances.items() if d == min_distance]
        candidate = ",".join(best_alleles)
        insertion_idx = None
        
        if candidate and candidate not in orig_call:
            merged, insertion_idx = merge_assignment_with_index(
                orig_call, candidate, sample_seq, germline, max_length=max_call
            )
            record[call] = merged
            
        # If a likelihood column is provided, update it.
        if likelihood_col and likelihood_col in record and insertion_idx is not None:
            curr_like = record[likelihood_col]
            
            # Remove brackets and split by whitespace
            values = curr_like.strip('[]').split()
            
            # Insert None at the appropriate position
            values.insert(insertion_idx, "None")
            
            # Trim if needed
            if len(values) > max_call:
                values = values[:max_call]
            
            # Reconstruct the string with brackets
            record[likelihood_col] = '[' + ' '.join(values) + ']'
            
    return records


def write_tsv(records, template, outfile, extra_fields):
    """Write records to a TSV file using AIRR's derive_rearrangement."""
    if not records:
        return
    writer = derive_rearrangement(outfile, template, fields=extra_fields)
    for r in records:
        writer.write(r)


# -- Main Function --

def reassign_alleles(input_path, out_name, germline_file, novel_file, 
                    call_column="v_call", seq_column="sequence_alignment",
                    likelihood_column=None, max_call=3):
    """
    Main function to reassign novel alleles in AIRR-formatted rearrangement records.
    
    Args:
        input_path: Path to input TSV file with rearrangement records
        out_name: Output file name prefix
        germline_file: Path to germline FASTA file
        novel_file: Path to novel alleles TSV file
        call_column: Column name for allele calls (default: v_call)
        seq_column: Column name for the aligned sequence (default: sequence_alignment)
        likelihood_column: Column name for likelihood values (optional)
        max_call: Maximum number of allele calls allowed per record (default: 3)
    """
    # Read input records using AIRR's read_rearrangement.
    print("Reading input records from:", input_path)
    reader = read_rearrangement(input_path)
    external_fields = reader.external_fields
    records = list(reader)
    
    # Read novel file using pandas.
    novel_df = pd.read_csv(novel_file, sep="\t", dtype=str)
    novel_alleles = novel_df['polymorphism_call'].dropna().unique().tolist()
    novel_genes = novel_df['gene'].dropna().unique().tolist()
    
    # Read germline FASTA file using Biopython.
    germline_records = list(SeqIO.parse(germline_file, "fasta"))
    germline_dict = {record.id: str(record.seq) for record in germline_records}
    
    # Build a regex from novel_genes to identify records of interest.
    novel_genes_regex = "|".join(novel_genes)
    
    # Identify records whose call (call_column) contains any novel gene.
    subset_indices = [i for i, r in enumerate(records) 
                     if re.search(novel_genes_regex, r.get(call_column, ""))]
    print(f"Found {len(subset_indices)} records with novel genes in column '{call_column}'.")
    
    # Extract all allele assignments from these records.
    novel_assignments_set = set()
    for i in subset_indices:
        call_str = records[i].get(call_column, "")
        for allele in call_str.split(","):
            novel_assignments_set.add(allele.strip())
    novel_assignments = list(novel_assignments_set)
    
    # Build germline_novel as the union of novel alleles (from the novel file) and assignments
    combined_alleles = set(novel_alleles) | set(novel_assignments)
    germline_novel = {allele: germline_dict[allele] 
                     for allele in combined_alleles if allele in germline_dict}
    print(f"Using {len(germline_novel)} alleles from combined novel data in germline_novel.")
    
    # Build a subset of records to update.
    records_to_update = [records[i] for i in subset_indices]
    updated_subset = reassign_novel_allele(
        records_to_update, 
        germline_novel,
        call=call_column,
        seq=seq_column,
        max_call=max_call,
        likelihood_col=likelihood_column
    )
    
    # Update the original records with the updated subset.
    for pos, idx in enumerate(subset_indices):
        records[idx] = updated_subset[pos]
    
    # Write output using derive_rearrangement.
    output_file = f"{out_name}_{call_column}_reassigned-pass.tsv"
    write_tsv(records, input_path, output_file, external_fields)
    print("Output written to", output_file)


# -- CLI --

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Novel Allele Reassignment (AIRR-compatible I/O)")
    parser.add_argument("--input", type=str, required=True,
                        help="Input TSV file with rearrangement records")
    parser.add_argument("--out_name", type=str, required=True,
                        help="Output file name prefix")
    parser.add_argument("--germline_file", type=str, required=True,
                        help="Germline FASTA file")
    parser.add_argument("--novel_file", type=str, required=True,
                        help="Novel alleles TSV file")
    parser.add_argument("--call_column", type=str, default="v_call",
                        help="Column name for allele calls (default: v_call)")
    parser.add_argument("--seq_column", type=str, default="sequence_alignment",
                        help="Column name for the aligned sequence (default: sequence_alignment)")
    parser.add_argument("--likelihood_column", type=str, default=None,
                        help="Column name for likelihood values (optional)")
    parser.add_argument("--max_call", type=int, default=3,
                        help="Maximum number of allele calls allowed per record (default: 3)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    reassign_alleles(
        input_path=args.input,
        out_name=args.out_name,
        germline_file=args.germline_file,
        novel_file=args.novel_file,
        call_column=args.call_column,
        seq_column=args.seq_column,
        likelihood_column=args.likelihood_column,
        max_call=args.max_call
    )
