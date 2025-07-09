#!/usr/bin/env python3

import argparse
import time
import itertools
from airr import read_rearrangement, derive_rearrangement

CHUNK_SIZE = 100_000

def clean_key(seq):
    """Remove gaps and periods from sequence_alignment."""
    return seq.replace('-', '').replace('.', '')

def pass1_minimal_collapse(args):
    """
    Pass 1: 
      - read the AIRR file in chunks
      - keep minimal data
      - track best row's 'sequence_id' for each KEY
      - sum consensus_count / duplicate_count
    Returns a dict collapsed_info: 
        KEY -> {
          'sequence_id': str,
          'consensus_count': int,
          'duplicate_count': int
        }
    """
    n_max = args.n_max
    reader = read_rearrangement(args.input)

    collapsed_info = {}  
    sums = {}            # sums[KEY] = {'consensus_count': int, 'duplicate_count': int}

    total_rows = 0
    chunk_index = 0

    while True:
        batch = list(itertools.islice(reader, CHUNK_SIZE))
        if not batch:
            break
        chunk_index += 1

        for row in batch:
            total_rows += 1
            vcall = row.get('v_call', '')
            jcall = row.get('j_call', '')
            if not (vcall and jcall):
                # skip if no v_call or j_call
                continue

            key = clean_key(row['sequence_alignment'])
            if key.count('N') > n_max:
                # skip high-N
                continue

            seq_id = row['sequence_id']
            cons = int(row.get('consensus_count', 0))
            dup  = int(row.get('duplicate_count', 0))

            if key not in collapsed_info:
                # store new best record
                collapsed_info[key] = {
                    'sequence_id': seq_id,
                    'n_count': key.count('N'),
                    'consensus_count': cons,
                    'duplicate_count': dup
                }
                sums[key] = {
                    'consensus_count': cons,
                    'duplicate_count': dup
                }
            else:
                # compare to existing best
                existing = collapsed_info[key]
                existing_n = existing['n_count']
                existing_cons = existing['consensus_count']

                # is this new row better?
                is_better = (key.count('N') < existing_n or 
                             cons > existing_cons)

                # always sum totals
                sums[key]['consensus_count'] += cons
                sums[key]['duplicate_count']  += dup

                if is_better:
                    # replace best row info
                    collapsed_info[key]['sequence_id'] = seq_id
                    collapsed_info[key]['n_count'] = key.count('N')
                    collapsed_info[key]['consensus_count'] = cons
                    collapsed_info[key]['duplicate_count'] = dup

    # inject final sums
    for key, info in collapsed_info.items():
        info['consensus_count'] = sums[key]['consensus_count']
        info['duplicate_count'] = sums[key]['duplicate_count']

    return collapsed_info


def pass2_write_full(args, collapsed_info):
    """
    Pass 2:
      - Re-read the input. 
      - For each row, see if row['sequence_id'] matches the best record for that KEY.
      - If it does, patch final sums & run pass/fail logic; write to outputs with full fields.
    """
    reader = read_rearrangement(args.input)
    conscount_flag = 'consensus_count' in reader.fields

    writer_pass = derive_rearrangement(f"{args.output_prefix}passed.tsv", args.input)
    writer_fail = derive_rearrangement(f"{args.output_prefix}failed.tsv", args.input)
    fasta_out_path = f"{args.output_prefix}seq.fasta"

    non_prod_removed = 0
    low_count_removed = 0
    pass_count = 0
    fail_count = 0
    ns_removed = 0  # Not strictly tracked here, done in pass1

    with open(fasta_out_path, 'w') as fasta_out:
        for row in reader:
            key = clean_key(row['sequence_alignment'])
            if key not in collapsed_info:
                # either not a valid key or filtered out
                continue

            # check if row matches the best row for that key
            best_info = collapsed_info[key]
            if row['sequence_id'] != best_info['sequence_id']:
                # not the chosen best
                continue

            # at this point, we have the row that is best
            # patch final sums
            row['consensus_count'] = best_info['consensus_count']
            row['duplicate_count']  = best_info['duplicate_count']

            # pass/fail logic
            if str(row.get('productive','')).upper() not in ['T','TRUE'] or 'J' not in row.get('j_call',''):
                non_prod_removed += 1
                writer_fail.write(row)
                fail_count += 1
            elif row.get('consensus_count' if conscount_flag else 'duplicate_count', 0) < args.min_conscount:
                low_count_removed += 1
                writer_fail.write(row)
                fail_count += 1
            else:
                pass_count += 1
                writer_pass.write(row)

                # create FASTA record
                # use KEY as 'sequence' or the original 'sequence' if you prefer
                seq = key  
                header_parts = []
                # if you want to add some fields, e.g. sequence_id, c_call, etc:
                for col in ['sequence_id', 'c_call', 'consensus_count', 'duplicate_count']:
                    if col in row and row[col]:
                        header_parts.append(f"{col}={row[col]}")
                header = '|'.join(header_parts).replace('sequence_id=', '')
                fasta_out.write(f">{header}\n{seq}\n")

    writer_pass.close()
    writer_fail.close()

    # summary
    return {
        'pass_count': pass_count,
        'fail_count': fail_count,
        'non_prod_removed': non_prod_removed,
        'low_count_removed': low_count_removed,
        'ns_removed': ns_removed
    }


def main():
    parser = argparse.ArgumentParser(description='Two-pass minimal-memory IGH collapse')
    parser.add_argument('--input', required=True, help='Input AIRR TSV file')
    parser.add_argument('--output_prefix', required=True, help='Output prefix')
    parser.add_argument('--n_max', type=int, default=10, help='Max allowed Ns before filtering out')
    parser.add_argument('--min_conscount', type=int, default=2, help='Min consensus_count threshold')
    parser.add_argument('--log_file', help='Optional log file')
    args = parser.parse_args()

    start = time.time()
    print('PASS1> Collapsing minimal data...')
    collapsed_info = pass1_minimal_collapse(args)
    print(f'PASS1> Found {len(collapsed_info)} unique alignment keys.\n')

    print('PASS2> Re-reading full file and writing final outputs...')
    stats = pass2_write_full(args, collapsed_info)
    duration = round((time.time() - start) / 60, 2)

    print('\nSUMMARY:')
    print(f"  Passed: {stats['pass_count']}")
    print(f"  Failed: {stats['fail_count']}")
    print(f"    - Non-productive or missing J: {stats['non_prod_removed']}")
    print(f"    - Below count threshold: {stats['low_count_removed']}")
    print("  (N removed was counted in pass1, not displayed here)")

    print(f"  Total time: {duration:.2f} min")
    print('DONE> Two-pass minimal memory collapse')

    if args.log_file:
        with open(args.log_file, 'w') as lf:
            lf.write(f"Passed: {stats['pass_count']}\n")
            lf.write(f"Failed: {stats['fail_count']}\n")
            lf.write(f"Non-productive or missing J: {stats['non_prod_removed']}\n")
            lf.write(f"Below count threshold: {stats['low_count_removed']}\n")
            lf.write(f"Unique keys: {len(collapsed_info)}\n")
            lf.write(f"Time (min): {duration}\n")


if __name__ == '__main__':
    main()
