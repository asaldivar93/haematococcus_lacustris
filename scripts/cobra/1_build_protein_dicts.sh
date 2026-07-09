#!/bin/bash
CWD=$HOME"/projects/h_lacustris"
GENOMES=$CWD"/data/external/genomes/h_lacustris"

RED_DB=$CWD"/data/external/genomes/h_lacustris/redball_algae_cb2023/redball.db"
HAEP_DB=$CWD"/data/external/genomes/h_lacustris/haeplu1/haeplu1.db"
NIES_DB=$CWD"/data/external/genomes/h_lacustris/nies_144/nies_144.db"

FMT_OUT="query,target,pident,evalue,bits,qcov,tcov,qlen,tlen"
FMT_CDD="query,target,pident,evalue,bits,qlen,tlen,qstart,qend,tstart,tend"

CDD_DB=$HOME"/databases/cdd/cdd.db"

RED_FA=$CWD"/data/external/genomes/h_lacustris/GCA_030144725.1_ASM3014472v1_protein.faa"
HAEP_FA=$CWD"/data/external/genomes/h_lacustris/GCA_050941845.1_Haeplu1_protein.faa"
NIES_FA=$CWD"/data/external/genomes/h_lacustris/GCA_011766145.1_Lacustris_1.0_protein.faa"

OUT_DIR=$CWD"/data/1_interim/genomes"

cd $CWD
mkdir -p $OUT_DIR

# mmseqs createdb haeplu1/*.faa haeplu1/haeplu1.db
# mmseqs createdb nies_144/*.faa nies_144/nies_144.db
# mmseqs createdb redball_algae_cb2023/*.faa redball_algae_cb2023/redball.db



# haep_to_red
mmseqs search $RED_DB $HAEP_DB \
	$GENOMES"/redball_algae_cb2023/haep_red.db" tmp \
	--threads 28 -s 7.5 -a
mmseqs convertalis $RED_DB $HAEP_DB \
	$GENOMES"/redball_algae_cb2023/haep_red.db" $OUT_DIR"haep_to_red_dict.tsv" \
	--format-output $FMT_OUT

# red to haep
mmseqs search $HAEP_DB $RED_DB \
	$GENOMES"/haeplu1/red_haep.db" tmp \
	--threads 28 -s 7.5 -a
mmseqs convertalis $RED_DB $HAEP_DB \
	$GENOMES"/haeplu1/red_haep.db" $OUT_DIR"red_to_haep_dict.tsv" \
	--format-output $FMT_OUT

# haep to nies
mmseqs search $NIES_DB $HAEP_DB \
	$GENOMES"/nies_144/haep_nies.db" tmp \
	--threads 28 -s 7.5 -a
mmseqs convertalis $NIES_DB $HAEP_DB \
	$GENOMES"/nies_144/haep_nies.db" $OUT_DIR"haep_to_nies_dict.tsv" \
	--format-output $FMT_OUT

# nies to haep
mmseqs search $HAEP_DB $NIES_DB \
	$GENOMES"/haeplu1/nies_haep.db" tmp \
	--threads 28 -s 7.5 -a
mmseqs convertalis $HAEP_DB $NIES_DB \
	$GENOMES"/haeplu1/nies_haep.db" $OUT_DIR"nies_to_haep_dict.tsv" \
	--format-output $FMT_OUT

# red to nies
mmseqs search $NIES_DB $RED_DB \
	$GENOMES"/nies_144/red_nies.db" tmp \
	--threads 28 -s 7.5 -a
mmseqs convertalis $NIES_DB $RED_DB \
	$GENOMES"/nies_144/red_nies.db" $OUT_DIR"red_to_nies_dict.tsv" \
	--format-output $FMT_OUT

#nies to red
mmseqs search $RED_DB $NIES_DB \
	$GENOMES"/redball_algae_cb2023/nies_red.db" tmp \
	--threads 28 -s 7.5 -a
mmseqs convertalis $RED_DB $NIES_DB \
	$GENOMES"/redball_algae_cb2023/nies_red.db" $OUT_DIR"nies_to_red_dict.tsv" \
	--format-output $FMT_OUT

#Reciprocal best hit
# mmseqs easy-rbh $RED_FA $HAEP_FA \
	# $OUT_DIR"/rbh_red_haep.tsv" tmp \
	# --threads 28 -s 7.5 \
	# --format-output $FMT_OUT

# mmseqs easy-rbh $RED_FA $NIES_FA \
	# $OUT_DIR"/rbh_red_nies.tsv" tmp \
	# --threads 28 -s 7.5 \
	# --format-output $FMT_OUT

# mmseqs easy-rbh $HAEP_FA $NIES_FA \
	# $OUT_DIR"/rbh_haep_nies.tsv" tmp \
	# --threads 28 -s 7.5 \
	# --format-output $FMT_OUT

#Annotate genomes with ncbi cdd
# mmseqs easy-search $RED_FA $CDD_DB \
	# $OUT_DIR"/cdd_red_hits.tsv" tmp \
	# --threads 28 -s 7.5 \
	# --format-output $FMT_CDD

# mmseqs easy-search $HAEP_FA $CDD_DB \
	# $OUT_DIR"/cdd_haep_hits.tsv" tmp \
	# --threads 28 -s 7.5 \
	# --format-output $FMT_CDD

# mmseqs easy-search $NIES_FA $CDD_DB \
	# $OUT_DIR"/cdd_nies_hits.tsv" tmp \
	# --threads 28 -s 7.5 \
	# --format-output $FMT_CDD

mmseqs easy-rbh $NIES_FA $RED_FA \
	$OUT_DIR"/rbh_nies_red.tsv" tmp \
	--threads 14 -s 7.5 \
	--format-output $FMT_OUT

mmseqs easy-rbh $NIES_FA $HAEP_FA \
	$OUT_DIR"/rbh_nies_haep.tsv" tmp \
	--threads 12 -s 7.5 \
	--format-output $FMT_OUT
