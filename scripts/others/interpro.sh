#!/bin/bash
awk '/^>/{sub(/[ \t].*/, ""); print; next} {print}' GCA_011766*.faa > nies.faa


nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile docker \
  --input data/external/genomes/h_lacustris/nies.faa \
  --datadir ~/databases/interpro \
  --interpro latest \
  --outdir data/2_processed/annotation/interpro \
  --outprefix nies144 \
  --skip-applications SMART \
  --goterms --pathways \
  --cpus 12 --max-workers 12

cut -f 1,2 nies144.tsv > nies144_md5_map.tsv

awk '/^>/{sub(/[ \t].*/, ""); print; next} {print}' GCA_030144725*.faa > ref.faa

nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile docker \
  --input data/external/genomes/h_lacustris/ref.faa \
  --datadir ~/databases/interpro \
  --interpro latest \
  --outdir data/2_processed/annotation/interpro \
  --outprefix ref \
  --skip-applications SMART \
  --goterms --pathways \
  --cpus 12 --max-workers 12
  
cut -f 1,2 ref.tsv > ref_md5_map.tsv
