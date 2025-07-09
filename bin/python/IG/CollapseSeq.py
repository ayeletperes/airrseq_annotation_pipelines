#!/usr/bin/env python3

"""
Collapse AIRR rearrangement records based on sequence and selected fields.

This tool is adapted from the Presto `CollapseSeq.py` logic but designed for
AIRR TSV files using the AIRR Python library.

https://github.com/immcantation/presto/blob/master/bin/CollapseSeq.py

Output files:
    - <output>_collapsed-pass.tsv: collapsed unique sequences.
    - <output>_collapsed-duplicate.tsv: duplicate sequences not retained.
    - <output>_collapsed-undetermined.tsv: ambiguous sequences over threshold.
    - <output>_collapsed-failed.tsv: failed initial filters.
    - <output>_collapsed-seq.fasta: optional FASTA output for downstream use.
    - <log_file>: log file with summary.

Output annotations:
    - duplicate_ids: IDs of all duplicate sequences collapsed.
    - User-defined fields (e.g., duplicate_count, consensus_count).
"""

import os
import re
import argparse
from collections import OrderedDict
from itertools import chain
from time import time
from airr import read_rearrangement, derive_rearrangement
from presto.Sequence import checkSeqEqual
from presto.IO import printProgress


# -- Constants and Defaults --
GAPP_RE = re.compile(r'[.]+')
AMBIG_RE = re.compile(r'[\-N]')
DEFAULT_MAX_MISSING = 10
LIKELIHOOD_RE = re.compile(r'\s+', re.I)


# -- Classes and Utilities --

class DuplicateSet:
    """Holds a collapsed sequence and its duplicates."""
    def __init__(self, row, key, missing, annotations):
        self.row = row
        self.missing = missing
        self.annotations = annotations
        self.keys = [key]
        self.count = 1

    def __len__(self):
        return len(self.keys)


def clean_gapps(seq):
    return GAPP_RE.sub('', seq)


def fill_missing_fields(row, n_max):
    """Sanitize a row and skip if exceeds ambiguity threshold."""
    if not (row.get('v_call') and row.get('j_call')):
        return None

    seq_alignment = row.get('sequence_alignment', '')
    seq_clean = clean_gapps(seq_alignment)
    if len(AMBIG_RE.findall(seq_clean)) > n_max:
        return None
    row['sequence_vdj'] = seq_clean

    if not row.get('c_call'):
        for alt in ['isotype', 'primer', 'reverse_primer', 'prcons', 'barcode']:
            if row.get(alt):
                row['c_call'] = row[alt]
                break
        else:
            row['c_call'] = 'Unk'

    row.setdefault('consensus_count', 1)
    row.setdefault('duplicate_count', 1)
    return row


def findUID(uid, search_dict, score=False):
    if not score:
        return uid if uid in search_dict else None
    for key in search_dict:
        if uid[1:] == key[1:] and checkSeqEqual(uid[0], key[0]):
            return key
    return None

def extract_first_likelihood(likelihood_str):
    """Extract the first value from a likelihood string like '[0.9244691 0.09451581]'"""
    if not likelihood_str:
        return 0.0
    try:
        # Remove brackets and split by whitespace
        values = likelihood_str.strip('[]').split()
        if not values:
            return 0.0
        return float(values[0])
    except (ValueError, IndexError):
        return 0.0

def find_unique_seq(uniq_dicts, search_keys, seq_dict, max_missing, uniq_fields, copy_fields, max_field=None, min_field=None):
    """Collapse duplicates by VDJ and selected fields."""
    dup_keys, to_remove = [], []
    start_time = time()

    for idx, key in enumerate(search_keys):
        printProgress(idx, len(search_keys), 0.05, start_time, task=f'{max_missing} missing')

        row = seq_dict[key]
        seq_str = row.get('sequence_vdj', '')
        ambig_count = len(AMBIG_RE.findall(seq_str))

        if ambig_count > max_missing:
            continue

        uid = tuple(chain([seq_str], (row.get(f, '') for f in uniq_fields)))
        cid = {f: [row.get(f, 0)] for f in copy_fields}
        seq_len = len(uid[0])
        uniq_dict = uniq_dicts.setdefault(seq_len, {})
        score = max_missing > 0

        match = findUID(uid, uniq_dict, score)
        if match is None:
            uniq_dict[uid] = DuplicateSet(row, key=key, missing=ambig_count, annotations=cid)
        else:
            entry = uniq_dict[match]
            entry.count += 1
            entry.keys.append(key)
            for f, v in cid.items():
                entry.annotations[f].extend(v)
            swap = False
            if ambig_count <= entry.missing:
                if max_field and row.get(max_field) and entry.row.get(max_field):
                    if 'likelihoods' in max_field.lower():
                            row_val = extract_first_likelihood(row[max_field])
                            entry_val = extract_first_likelihood(entry.row[max_field])
                            swap = row_val > entry_val
                    else:
                        swap = float(row[max_field]) > float(entry.row[max_field])
                elif min_field and row.get(min_field) and entry.row.get(min_field):
                    # Similar handling for min_field if needed
                    if 'likelihoods' in min_field.lower():
                            row_val = extract_first_likelihood(row[min_field])
                            entry_val = extract_first_likelihood(entry.row[min_field])
                            swap = row_val < entry_val
                    else:
                            swap = float(row[min_field]) < float(entry.row[min_field])
            if swap:
                dup_keys.append(entry.row['sequence_id'])
                entry.row = row
                entry.missing = ambig_count
            else:
                dup_keys.append(key)
        to_remove.append(idx)
    for j in reversed(to_remove):
        del search_keys[j]
    printProgress(len(search_keys), len(search_keys), 0.05, start_time, task=f'{max_missing} missing')
    return uniq_dicts, search_keys, dup_keys


def merge_dicts(dicts):
    return {k: v for _, d in dicts.items() for k, v in d.items()}

def collapse_likelihoods(likelihoods):
    """
    Process likelihood strings by converting space-separated values to comma-separated.
    
    Args:
        likelihoods: List of likelihood strings like ['[0.9244691 0.09451581]']
    
    Returns:
        String representation of the values joined with commas
    """
    result = []
    for likelihood_str in likelihoods:
        # Remove brackets and split by whitespace
        values = likelihood_str.strip('[]').split()
        # Join with commas
        comma_joined = ','.join(values)
        result.append(comma_joined)
    
    return ','.join(result)

def collapse_airr_annotations(ann_dict, fields, actions, delimiter=','):
    """
    Collapse annotations based on specified actions.
    
    Args:
        ann_dict: Dictionary of annotations
        fields: List of field names
        actions: List of actions to apply
        delimiter: Delimiter for 'set' action
    
    Returns:
        Dictionary with collapsed annotations
    """
    def apply(values, action):
        if not values:
            return None
        
        if action == 'sum':
            return str(sum(float(x or 0) for x in values))
        elif action == 'sum_int':
            return sum(int(x or 0) for x in values)
        elif action == 'min':
            return f"{min(float(x or 0) for x in values):.12g}"
        elif action == 'max':
            return f"{max(float(x or 0) for x in values):.12g}"
        elif action == 'first':
            return values[0]
        elif action == 'last':
            return values[-1]
        elif action == 'set':
            return delimiter.join(sorted(set(str(x) for x in values if x)))
        elif action == 'cat':
            return ''.join(str(x) for x in values)
        elif action == 'likelihood':
            return collapse_likelihoods(values)
        
        return values
    
    return {f: apply(ann_dict.get(f, []), a) for f, a in zip(fields, actions)}

def write_airr_records(records, output_path, input_path, include_fields):
    if not records:
        return
    writer = derive_rearrangement(output_path, input_path, fields=include_fields)
    for row in records:
        writer.write(row)


def write_log_file(input_path, output_path, uniq_dict, dup_keys, search_keys, seq_dict,
                   max_missing, uniq_fields, copy_fields, copy_actions):
    log = OrderedDict()
    log['START'] = 'CollapseSeq'
    log['FILE'] = os.path.basename(input_path)
    log['MAX_MISSING'] = max_missing
    log['UNIQ_FIELDS'] = ','.join(uniq_fields)
    log['COPY_FIELDS'] = ','.join(copy_fields)
    log['COPY_ACTIONS'] = ','.join(copy_actions)
    log['OUTPUT'] = os.path.basename(output_path)
    log['SEQUENCES'] = len(seq_dict)
    log['UNIQUE'] = len(uniq_dict)
    log['DUPLICATE'] = len(dup_keys)
    log['UNDETERMINED'] = len(search_keys)
    log['END'] = 'CollapseSeq'

    path = output_path.replace(".tsv", ".log")
    with open(path, "w") as f:
        for k, v in log.items():
            f.write(f"{k}: {v}\n")


# -- Main Collapse Function --

def collapse_seq(input_path, output_prefix, max_missing, uniq_fields, copy_fields, copy_actions,
                 max_field=None, min_field=None, remove_consensus_count=False,
                 write_failed=True,write_duplicates=True,log_file=None):

    reader = read_rearrangement(input_path)
    external_fields = reader.external_fields
    external_actions = []
    for field in external_fields:
        if 'likelihoods' in field.lower():
            continue
        else:
            external_actions.append('set')
            copy_fields.append(field)
            
    copy_actions = copy_actions + external_actions

    seq_dict, failed_records = {}, []
    for rec in reader:
        filled = fill_missing_fields(rec, max_missing)
        if filled:
            seq_dict[filled['sequence_id']] = filled
        else:
            failed_records.append(rec)

    search_keys = list(seq_dict.keys())
    uniq_dicts, dup_keys = {}, []

    if max_field:
        print(f"Using {max_field} field for choosing the collapse sequence.")
        
    for n in range(0, max_missing + 1):
        uniq_dicts, search_keys, dup_list = find_unique_seq(
            uniq_dicts, search_keys, seq_dict, n,
            uniq_fields=uniq_fields,
            copy_fields=copy_fields,
            max_field=max_field,
            min_field=min_field
        )
        dup_keys.extend(dup_list)
        if not search_keys:
            break

    uniq_dict = merge_dicts(uniq_dicts)

    # --- Write outputs ---
    collapsed_records, collapsed_fasta, duplicate_records = [], [], []
    for entry in uniq_dict.values():
        row = entry.row.copy()
        agg = collapse_airr_annotations(entry.annotations, copy_fields, copy_actions)
        row.update(agg)
        row['duplicate_ids'] = ','.join(entry.keys)
        collapsed_records.append(row)
        duplicate_records.extend(seq_dict[k] for k in entry.keys[1:])

        fasta_id = row['sequence_id'] + "|" + "|".join(
            f"{f}={row.get(f, '')}" for f in (copy_fields + ['c_call']) if f in row
        )
        collapsed_fasta.append({'id': fasta_id, 'seq': row.get('sequence', '')})

    if remove_consensus_count:
        for row in collapsed_records:
            row.pop('consensus_count', None)

    if collapsed_records:
        write_airr_records(collapsed_records, f"{output_prefix}_collapsed-passed.tsv", input_path,
                           include_fields=copy_fields + ['duplicate_ids'])

    if write_failed and failed_records:
        write_airr_records(failed_records, f"{output_prefix}_collapsed-failed.tsv", input_path,
                           include_fields=copy_fields)

    if write_duplicates and duplicate_records:
        write_airr_records(duplicate_records, f"{output_prefix}_collapsed-duplicate.tsv", input_path,
                           include_fields=copy_fields)

    if search_keys:
        undetermined = [seq_dict[k] for k in search_keys]
        write_airr_records(undetermined, f"{output_prefix}_collapsed-undetermined.tsv", input_path,
                           include_fields=copy_fields)

    if collapsed_fasta:
        with open(f"{output_prefix}_collapsed-seq.fasta", "w") as fh:
            for rec in collapsed_fasta:
                fh.write(f">{rec['id']}\n{rec['seq']}\n")

    if log_file:
        write_log_file(input_path, log_file, uniq_dict, dup_keys, search_keys, seq_dict, max_missing,
                       uniq_fields, copy_fields, copy_actions)


# -- CLI --

def parse_args():
    parser = argparse.ArgumentParser(description="Collapse AIRR rearrangement TSV files.")
    parser.add_argument("-i", "--input", required=True, help="Input AIRR TSV file.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix.")
    parser.add_argument("-n", "--max-missing", type=int, default=DEFAULT_MAX_MISSING,
                        help="Max number of N/- allowed per sequence.")
    parser.add_argument("--uf", nargs='+', default=['v_call', 'j_call', 'c_call'],
                        help="Fields that define uniqueness.")
    parser.add_argument("--cf", nargs='+', default=['duplicate_count', 'consensus_count'],
                        help="Fields to aggregate from duplicates.")
    parser.add_argument("--act", nargs='+', default=['sum_int', 'sum_int'],
                        help="Actions for copy fields (sum, set, min, max).")
    parser.add_argument("--maxf", dest="max_field", help="Field whose MAX value wins representative.")
    parser.add_argument("--minf", dest="min_field", help="Field whose MIN value wins representative.")
    parser.add_argument("--remove-consensus-count", action="store_true",
                   help="Drop consensus_count column from final TSV.")
    parser.add_argument("--no-failed",     action="store_false", dest="write_failed",
                   help="Do not write *_failed.tsv.")
    parser.add_argument("--no-duplicates", action="store_false", dest="write_duplicates",
                   help="Do not write *_duplicate.tsv.")
    parser.add_argument("--log-file", type=str, default=None,
                        help="Log file to write summary of filtering results")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    collapse_seq(
        input_path=args.input,
        output_prefix=args.output_prefix,
        max_missing=args.max_missing,
        uniq_fields=args.uf,
        copy_fields=args.cf,
        copy_actions=args.act,
        max_field=args.max_field,
        min_field=args.min_field,
        remove_consensus_count=args.remove_consensus_count,
        write_failed=args.write_failed,
        write_duplicates=args.write_duplicates,
        log_file=args.log_file
    )
