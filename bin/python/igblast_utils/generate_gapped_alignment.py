
#!/usr/bin/env python3

import argparse
import time
from airr import read_rearrangement, derive_rearrangement
from concurrent.futures import ProcessPoolExecutor
from itertools import islice, tee
import multiprocessing

# Global germline maps
germline_refs = {'v': {}, 'd': {}, 'j': {}}
IS_HEAVY = True

def parse_fasta(file):
    header, seq = None, []
    for line in file:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if header:
                yield header, ''.join(seq)
            header, seq = line, []
        else:
            seq.append(line)
    if header:
        yield header, ''.join(seq)

def load_germline_sequences(filepath):
    ref_map = {}
    with open(filepath) as f:
        for header, seq in parse_fasta(f):
            allele = header.lstrip('>').split('|')[0]
            ref_map[allele] = seq
    return ref_map

def parse_sequence_id(row):
    sid = row.get("sequence_id", "")
    parts = sid.split('|')
    row['sequence_id'] = parts[0]
    for p in parts[1:]:
        if '=' in p:
            k, v = p.split('=', 1)
            if k == 'CONSCOUNT':
                row['consensus_count'] = int(v)
            elif k == 'DUPCOUNT':
                row['duplicate_count'] = int(v)
            else:
                row[k.lower()] = v
    return row

def get_sequence_alignment(row):
    v_call = row.get('v_call', '').split(',')[0]
    print(v_call)
    v_ref = germline_refs['v'].get(v_call)
    print(v_ref)
    if not v_ref:
        return ''
    v_ref = v_ref[:row.get('v_germline_end', 0)]
    print(v_ref)
    seq = row['sequence'][row.get('v_sequence_start', 0):row.get('j_sequence_end', 0)]
    print(seq)
    print(row.get('v_germline_start', 0))
    if row.get('v_germline_start', 0) > 0:
        seq = '.' * row['v_germline_start'] + seq
    aligned = []
    seq_iter = iter(seq)
    started = False
    for base in v_ref:
        if base != '.':
            started = True
            aligned.append(next(seq_iter, '.'))
        else:
            aligned.append(next(seq_iter, '.') if started else '.')
    print(aligned)
    aligned.extend(seq_iter)
    return ''.join(aligned)

def get_germline_alignment(row):
    v_call = row.get('v_call', '').split(',')[0]
    j_call = row.get('j_call', '').split(',')[0]
    v_ref = germline_refs['v'].get(v_call, '')[:row.get('v_germline_end', 0)]
    j_ref = germline_refs['j'].get(j_call, '')[row.get('j_germline_start', 0):row.get('j_germline_end', 0)]
    d_ref = ''
    if IS_HEAVY:
        d_call = row.get('d_call', '').split(',')[0]
        d_raw = germline_refs['d'].get(d_call, '')
        d_ref = (
            row['sequence'][row['v_sequence_end']:row['d_sequence_start']] +
            d_raw[row.get('d_germline_start', 0):row.get('d_germline_end', 0)] +
            row['sequence'][row['d_sequence_end']:row['j_sequence_start']]
        )
    else:
        d_ref = row['sequence'][row['v_sequence_end']:row['j_sequence_start']]
    return v_ref + d_ref + j_ref

def process_chunk(chunk):
    output_rows = []
    for row in chunk:
        row = parse_sequence_id(row)
        row['sequence_alignment'] = get_sequence_alignment(row)
        row['germline_alignment'] = get_germline_alignment(row)
        output_rows.append(row)
    return output_rows

def chunked_iterator(iterable, size):
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk

def main():
    global germline_refs, IS_HEAVY

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--v_germline', required=True)
    parser.add_argument('--d_germline', help='D germline FASTA (for heavy chains)')
    parser.add_argument('--j_germline', required=True)
    parser.add_argument('--is_light', action='store_true')
    parser.add_argument('--threads', type=int, default=multiprocessing.cpu_count())
    parser.add_argument('--chunk_size', type=int, default=5000)
    parser.add_argument('--log_file', help='Optional log file')
    args = parser.parse_args()

    start = time.time()
    IS_HEAVY = not args.is_light
    germline_refs['v'] = load_germline_sequences(args.v_germline)
    germline_refs['j'] = load_germline_sequences(args.j_germline)
    if args.d_germline and IS_HEAVY:
        germline_refs['d'] = load_germline_sequences(args.d_germline)

    reader = read_rearrangement(args.input)
    duration = round(time.time() - start, 2)
    print(f"Read rows in {duration} seconds.")
    original_fields = set(reader.fields)
    
    all_rows = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(process_chunk, chunked_iterator(reader, args.chunk_size))
        for chunk in results:
            all_rows.extend(chunk)

    duration = round(time.time() - start, 2)
    print(f"Processed chunks in {duration} seconds.")
    
    all_output_keys = set()
    for row in all_rows:
        all_output_keys.update(row.keys())
    new_fields = sorted([f for f in (all_output_keys - original_fields) if f != 'KEY'])
    print(f"New fields: {new_fields}")
    
    writer = derive_rearrangement(args.output, args.input, fields=new_fields)
    for row in all_rows:
        writer.write(row)

    duration = round(time.time() - start, 2)
    print(f"Processed {len(all_rows)} rows in {duration} seconds.")
    if args.log_file:
        with open(args.log_file, 'w') as log:
            log.write(f"Processed rows: {len(all_rows)}\n")
            log.write(f"Elapsed time: {duration} seconds\n")

if __name__ == '__main__':
    main()
