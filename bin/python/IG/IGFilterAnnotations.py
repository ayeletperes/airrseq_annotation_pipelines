#!/usr/bin/env python3
"""
IG filter Seq: filter airr-seq annotation before downstream analysis.

Filters:
  - Remove sequences with excessive N in the V‑region.
  - Remove non-functional sequences (not productive, not VJ in-frame, or with stop codon).
  - Remove sequences with junction_length ≤ 0.
  - Remove sequences with low consensus_count (if field present).

Outputs:
  <prefix>_filtered-passed.tsv                      -- sequences passing all filters
  <prefix>_filtered-failed.tsv                      -- failed productive/in-frame checks; failed junction_length check; failed J call check; failed min consensus_count filter
  <prefix>_filtered-seq.fasta                       -- optional FASTA output
  <log_file>                                        -- log file with summary.
"""

import argparse
import re
from airr import read_rearrangement, derive_rearrangement
from presto.IO import printProgress
from time import time

N_RE = re.compile(r'N', re.I)
DEFAULT_MAX_MISSING = 10
DEFAULT_MIN_CONSENSUS = 0

def v_region_ok(row, n_limit):
    """
    Return True if number of Ns in V‑region (up to v_germline_end) is within threshold.
    """
    v_end = int(row.get("v_germline_end", 0))
    if v_end <= 1:
        return True
    return N_RE.findall(row["sequence_alignment"][: v_end - 1]).count("N") <= n_limit

def write_tsv(records, template, outfile, extra_fields):
    if not records:
        return
    writer = derive_rearrangement(outfile, template, fields=extra_fields)
    for r in records:
        writer.write(r)

def write_fasta(records, header_fields, outfile):
    with open(outfile, "w") as fh:
        for r in records:
            fasta_id = r['sequence_id'] + "|" + "|".join(
                f"{f}={r.get(f, '')}" for f in (header_fields) if f in r
            )
            fh.write(f">{fasta_id}\n{r.get('sequence','')}\n")

def filter_annotations(args):
    
    annotations_count = sum(1 for _ in read_rearrangement(args.input))
    reader = read_rearrangement(args.input)
    external_fields = reader.external_fields
    copy_fields = ['consensus_count', 'duplicate_count', 'c_call'] + external_fields
    
    
    seq_seen = set()
    passed = []
    non_func = []
    low_cons = []
    invalid_junction = []
    invalid_j_call = []
    dup_drop = []
    ns_removed = 0
    start_time = time()
    for idx, row in enumerate(reader):
        printProgress(idx, annotations_count, 0.05, start_time, task='Processing Annotations')
        
        sid = row["sequence_id"]
        if sid in seq_seen:
            dup_drop.append(row)
            continue
        seq_seen.add(sid)
        
        # V region N filter
        if not v_region_ok(row, args.n_limit):
            ns_removed+= 1
            continue

        # Productive filtering
        prod = row.get("productive", "")
        in_frame = row.get("vj_in_frame", "")
        stop_codon = row.get("stop_codon", "")

        if not prod or not in_frame or stop_codon:
            non_func.append(row)
            continue

        # Junction length check
        if int(row.get("junction_length", 0)) <= 0:
            invalid_junction.append(row)
            continue

        # J call check
        if not row.get("j_call", ""):
            invalid_j_call.append(row)
            continue
        
        # Consensus count filter
        if "consensus_count" in row and args.min_consensus > 0:
            try:
                if float(row["consensus_count"]) < args.min_consensus:
                    low_cons.append(row)
                    continue
            except ValueError:
                low_cons.append(row)
                continue

        passed.append(row)

    printProgress(annotations_count, annotations_count, 0.05, start_time, task='Processing Annotations')
    # combine non valid sequences to a single list
    non_func.extend(invalid_junction)
    non_func.extend(invalid_j_call)
    non_func.extend(dup_drop)
    non_func.extend(low_cons)
    
    prefix = args.output_prefix
    write_tsv(passed, args.input, f"{prefix}_filtered-passed.tsv", external_fields)
    

    if not args.write_failed:
        write_tsv(non_func, args.input, f"{prefix}_filtered-failed.tsv", external_fields)
            
    if args.fasta:
        write_fasta(passed, copy_fields, f"{prefix}_filtered-seq.fasta")
    
    if args.log_file:
        with open(args.log_file, 'w') as log:
            log.write(f"NS removed: {ns_removed}\n")
            log.write(f"Non-productive: {len(non_func)}\n")
            log.write(f"Missing J removed: {len(invalid_j_call)}\n")
            log.write(f"Invalid junction: {len(invalid_junction)}\n")
            if args.min_consensus > 0:
                log.write(f"Low consensus removed: {len(low_cons)}\n")
            log.write(f"Final sequences: {len(passed)}\n")
            
def parse_args():
    parser = argparse.ArgumentParser(description="Filter IG annotation table prior to collapse.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output_prefix", required=True)
    parser.add_argument("--n-limit", type=int, default=DEFAULT_MAX_MISSING,
                        help="Max number of Ns allowed in the V-region (default: 10)")
    parser.add_argument("--min-consensus", type=float, default=DEFAULT_MIN_CONSENSUS,
                        help="Minimum consensus_count required (default: 0)")
    parser.add_argument("--fasta", action="store_true",
                        help="Also write FASTA file from passed sequences")
    parser.add_argument("--log-file", type=str, default=None,
                        help="Log file to write summary of filtering results")
    parser.add_argument("--no-failed",     action="store_false", dest="write_failed",
                   help="Do not write *_failed.tsv.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    filter_annotations(args)
