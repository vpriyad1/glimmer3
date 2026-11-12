#!/usr/bin/env bash
# Run Glimmer3 on pO157 (AB011549.2) and score it against the GenBank annotation.
# Requires: glimmer3, long-orfs, extract, build-icm on PATH (conda install -c bioconda glimmer)
set -euo pipefail

ACC="AB011549.2"
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&retmode=text"
mkdir -p data glimmer_run results

# 1. Download the reference record and sequence (skipped if already present)
[ -s "data/${ACC}.gbk" ]   || wget -q -O "data/${ACC}.gbk"   "${EFETCH}&rettype=gbwithparts"
[ -s "data/${ACC}.fasta" ] || wget -q -O "data/${ACC}.fasta" "${EFETCH}&rettype=fasta"

# 2. Train Glimmer3 on the sequence's own long ORFs, then predict
long-orfs -n -t 1.15 "data/${ACC}.fasta" glimmer_run/longorfs.coords
extract -t "data/${ACC}.fasta" glimmer_run/longorfs.coords > glimmer_run/train.seq
build-icm -r glimmer_run/glimmer.icm < glimmer_run/train.seq
glimmer3 "data/${ACC}.fasta" glimmer_run/glimmer.icm "data/${ACC}_glimmer"

# 3. Parse both gene sets and compare them
SEQ_LENGTH=$(grep -v "^>" "data/${ACC}.fasta" | tr -d '\n' | wc -c)
python3 scripts/parse_gbk.py "data/${ACC}.gbk" results/reference_cds.json
python3 scripts/parse_glimmer.py "data/${ACC}_glimmer.predict" results/glimmer_predictions.json --length "${SEQ_LENGTH}"
python3 scripts/compare_ref_vs_glimmer.py results/reference_cds.json results/glimmer_predictions.json results
