#!/usr/bin/env bash

# Usage: ./reverse_rename_fasta_ids.sh changes.csv annotation.tsv germline.fasta

CHANGES_CSV="$1"
ANNOTATION_FILE="$2"
GERMLINE_FILE="$3"

if [ ! -f "$CHANGES_CSV" ]; then
  echo "No changes.csv file found. Skipping rename reversal."
  exit 0
fi

while IFS=, read -r _ OLD_ID NEW_ID; do
  OLD_SUFFIX="${OLD_ID#*\*}"
  NEW_SUFFIX="${NEW_ID#*\*}"

  # Replace in annotation and germline files
  sed -i "s/${NEW_SUFFIX}/${OLD_SUFFIX}/g" "$ANNOTATION_FILE"
  sed -i "s/${NEW_SUFFIX}/${OLD_SUFFIX}/g" "$GERMLINE_FILE"

done < <(tail -n +2 "$CHANGES_CSV")
