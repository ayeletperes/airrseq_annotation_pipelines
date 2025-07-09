#!/usr/bin/env python3

"""
Script to clean and rename long FASTA IDs using SHA256 suffix,
while preserving unique identifiers. Also outputs a CSV log
of any ID changes.
"""

from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
from hashlib import sha256
import sys
import csv

input_fasta, output_fasta, output_changes_csv = sys.argv[1:]

records = list(SeqIO.parse(input_fasta, "fasta"))
new_records = []
changes = []

for i, rec in enumerate(records, 1):
    old_id = rec.id
    if len(old_id) > 50 and '*' in old_id:
        gene, rest = old_id.split('*', 1)
        allele = rest.split('_')[0]
        suffix = sha256(str(rec.seq).encode()).hexdigest()[-4:]
        new_id = f"{gene}*{allele}_{suffix}"
        rec.id = new_id
        rec.description = ''
        if old_id != new_id:
            changes.append([i, old_id, new_id])
    new_records.append(rec)

SeqIO.write(new_records, output_fasta, "fasta")

with open(output_changes_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Row", "Old_ID", "New_ID"])
    writer.writerows(changes)
