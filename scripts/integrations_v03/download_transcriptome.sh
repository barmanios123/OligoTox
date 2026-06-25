#!/bin/bash
set -e
mkdir -p data/raw/transcriptome
echo "Downloading Ensembl GRCh38 cDNA FASTA (release 111, ~350MB)..."
wget -nc -q \
  "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
  -O data/raw/transcriptome/ensembl_human_cdna_111.fa.gz
echo "Done. Set use_real_transcriptome: true in configs/integration_config.yaml"
