#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import argparse
import time
from collections import Counter
import airr

SPLITSIZE = 2

class CollapseDict:
    def __init__(self, iterable=(), depth=0, nlim=10, conscount_flag=False):
        self.lowqual = {}
        self.seqs = {}
        self.children = {}
        self.depth = depth
        self.nlim = nlim
        self.conscount = conscount_flag
        for record in iterable:
            self.add(record)

    def split(self):
        newseqs = {}
        for seq in self.seqs:
            if len(seq) == self.depth:
                newseqs[seq] = self.seqs[seq]
            else:
                if seq[self.depth] not in self.children:
                    self.children[seq[self.depth]] = CollapseDict(depth=self.depth + 1)
                self.children[seq[self.depth]].add(self.seqs[seq], seq)
        self.seqs = newseqs

    def add(self, record, key=None):
        key = record['sequence_alignment'].replace('-', '').replace('.', '')
        record['KEY'] = key
        record.setdefault('ISOTYPECOUNTER', Counter([record.get('c_call', 'Ig')]))
        record.setdefault('VGENECOUNTER', Counter([record.get('v_call', '')]))
        record.setdefault('JGENECOUNTER', Counter([record.get('j_call', '')]))
        key = key or record['KEY']

        if self.depth == 0:
            if not record.get('j_call') or not record.get('v_call'):
                return
            if key.count('N') > self.nlim:
                self.lowqual[key] = combine(self.lowqual.get(key, {}), record, self.conscount)
                return

        if len(self.seqs) > SPLITSIZE:
            self.split()

        if key in self.seqs:
            self.seqs[key] = combine(self.seqs[key], record, self.conscount)
        elif key[self.depth] in self.children:
            self.children[key[self.depth]].add(record, key)
        else:
            self.seqs[key] = record

    def __iter__(self):
        yield from self.seqs.items()
        for child in self.children.values():
            yield from child
        yield from self.lowqual.items()

def combine(r1, r2, conscount_flag):
    def val(r): return -r['KEY'].count('N'), int(r.get('consensus_count', 0))
    if not r1:
        return r2
    result = r1 if val(r1) >= val(r2) else r2.copy()
    if conscount_flag:
        result['consensus_count'] = int(r1.get('consensus_count', 0)) + int(r2.get('consensus_count', 0))
    result['duplicate_count'] = int(r1.get('duplicate_count', 0)) + int(r2.get('duplicate_count', 0))
    result['ISOTYPECOUNTER'] = r1.get('ISOTYPECOUNTER', Counter()) + r2.get('ISOTYPECOUNTER', Counter())
    result['VGENECOUNTER'] = r1.get('VGENECOUNTER', Counter()) + r2.get('VGENECOUNTER', Counter())
    result['JGENECOUNTER'] = r1.get('JGENECOUNTER', Counter()) + r2.get('JGENECOUNTER', Counter())
    return result

def print_progress_bar(current, total, bar_length=24):
    percent = current / total
    blocks = int(bar_length * percent)
    bar = '#' * blocks + ' ' * (bar_length - blocks)
    print(f"PROGRESS> {time.strftime('%H:%M:%S')} |{bar}| {int(percent*100)}% ({current})", flush=True)

def main():
    parser = argparse.ArgumentParser(description="IGH AIRR-seq collapse and filtering script")
    parser.add_argument("--input", required=True, help="Input AIRR TSV file")
    parser.add_argument("--output_prefix", required=True, help="Output file prefix")
    parser.add_argument("--n_max", type=int, default=10, help="Maximum allowed Ns")
    parser.add_argument("--min_conscount", type=int, default=2, help="Minimum consensus count")
    parser.add_argument("--log_file", help="Optional log file for statistics")
    args = parser.parse_args()

    print(f"   START> IGCollapse")
    start_time = time.time()
    print(f"PROGRESS> {time.strftime('%H:%M:%S')} |Loading files       | 0.0 min")

    
    airrseq_data = pd.read_csv(args.input, sep='\t')
    airrseq_data['sequence_vdj'] = airrseq_data['sequence_alignment'].str.replace('[-.]', '', regex=True)
    airrseq_data['KEY'] = airrseq_data['sequence_vdj']
    airrseq_data['duplicate_count'] = airrseq_data.get('duplicate_count', 1)

    if 'c_call' not in airrseq_data.columns:
        for alt_col in ['isotype', 'primer', 'reverse_primer', 'prcons', 'barcode']:
            if alt_col in airrseq_data.columns:
                airrseq_data['c_call'] = airrseq_data[alt_col]
                break
        else:
            airrseq_data['c_call'] = 'Ig'

    airrseq_data.drop_duplicates(subset='sequence_id', inplace=True)
    pre_filter_ns = len(airrseq_data)
    v_end = airrseq_data['v_germline_end'].fillna(0).astype(int)
    airrseq_data['N_count'] = [
        seq[:ve].count('N') if isinstance(seq, str) else 0
        for seq, ve in zip(airrseq_data['sequence_alignment'], v_end)
    ]
    airrseq_data = airrseq_data[airrseq_data['N_count'] <= args.n_max]
    ns_removed = pre_filter_ns - len(airrseq_data)

    conscount_flag = 'consensus_count' in airrseq_data.columns
    records = airrseq_data.to_dict(orient='records')
    airrseq_data_collapsed = CollapseDict(records, nlim=args.n_max, conscount_flag=conscount_flag)

    collapsed_records = []
    for i, (_, rec) in enumerate(airrseq_data_collapsed):
        rec['sequence'] = rec['KEY']
        if conscount_flag:
            rec['consensus_count'] = int(rec.get('consensus_count', 0))
        rec['duplicate_count'] = int(rec.get('duplicate_count', 0))
        collapsed_records.append(rec)
        if i % max(1, len(records)//20) == 0:
            print_progress_bar(i, len(records))

    df_passed = pd.DataFrame(collapsed_records)
    df_passed.drop(columns=['sequence_vdj', 'KEY'], errors='ignore', inplace=True)

    pre_productive = len(df_passed)
    df_passed = df_passed[
        df_passed['productive'].astype(str).str.upper().isin(['T', 'TRUE']) &
        df_passed['j_call'].astype(str).str.contains('J')
    ]
    non_prod_removed = pre_productive - len(df_passed)

    filter_col = 'consensus_count' if conscount_flag else 'duplicate_count'
    df_low = df_passed[df_passed[filter_col] < args.min_conscount]
    df_passed = df_passed[df_passed[filter_col] >= args.min_conscount]

    df_passed.to_csv(f"{args.output_prefix}passed.tsv", sep='\t', index=False)
    
    #airr.dump_rearrangement(df_passed, f"{args.output_prefix}passed.tsv")
    
    pd.concat([df_low]).to_csv(f"{args.output_prefix}failed.tsv", sep='\t', index=False)

    #airr.dump_rearrangement(pd.concat([df_low]), f"{args.output_prefix}failed.tsv")
    
    # Write FASTA file using only relevant metadata columns (cdr3+1 onward)
    fasta_output_path = f"{args.output_prefix}seq.fasta"
    all_columns = list(df_passed.columns)
    core_fields = ['sequence_id', 'duplicate_count', 'consensus_count', 'c_call']
    if 'cdr3' in all_columns:
        idx_cdr3 = all_columns.index('cdr3') + 1
        post_cdr3_columns = all_columns[idx_cdr3:]
    else:
        post_cdr3_columns = []

    unique_information = list(dict.fromkeys(core_fields + post_cdr3_columns))
    unique_information = [col for col in unique_information if col in all_columns]

    with open(fasta_output_path, 'w') as fasta_out:
        for _, row in df_passed.iterrows():
            header_parts = []
            for col in unique_information:
                if pd.notna(row.get(col, None)):
                    header_parts.append(f"{col}={row[col]}")
            header = '|'.join(header_parts).replace('sequence_id=', '')
            fasta_out.write(f">{header}\n{row['sequence']}\n")

    total_time = round((time.time() - start_time) / 60, 2)
    print(f"PROGRESS> {time.strftime('%H:%M:%S')} |Done                | {total_time:.1f} min")
    print(f"  OUTPUT> {args.output_prefix}passed.tsv")
    print(f"    PASS> {len(df_passed)}")
    print(f"    FAIL> {len(df_low)}")
    print(f"     END> IGCollapse")

    if args.log_file:
        with open(args.log_file, 'w') as log:
            log.write(f"NS removed: {ns_removed}\n")
            log.write(f"Non-productive or missing J removed: {non_prod_removed}\n")
            log.write(f"Low {filter_col} removed: {len(df_low)}\n")
            log.write(f"Final sequences: {len(df_passed)}\n")

if __name__ == "__main__":
    main()
