import polars as pl

proteomics_lazy = pl.scan_csv(
    "data/3_results/proteomics_map.csv"
)
proteomics_lazy.with_columns(pl.col("Gene").str.split(";")).collect()

gene_info = pl.scan_csv(
    "data/external/genomes/h_lacustris/ncbi_dataset.tsv", separator="\t"
)

red_lazy = pl.scan_csv(
    "data/1_interim/genomes/rbh_nies_red.tsv",
    separator = "\t",
    has_header = False,
    new_columns = ["query","red","pident","evalue","bits","qcov","tcov","qlen","tlen"]
)

haep_lazy = pl.scan_csv(
    "data/1_interim/genomes/rbh_nies_haep.tsv",
    separator = "\t",
    has_header = False,
    new_columns = ["query","haep","pident","evalue","bits","qcov","tcov","qlen","tlen"]
)
haep_lazy.collect()

query = (
    proteomics_lazy
    .with_columns(pl.col("Gene").str.split(";"))
    .explode("Gene")
    .join(gene_info.select("Locus tag", "Protein accession"), left_on="Gene", right_on="Locus tag", how="left")
    .join(red_lazy.select("query", "red"), left_on="Protein accession", right_on="query", how="left")
    .join(haep_lazy.select("query", "haep"), left_on="Protein accession", right_on="query", how="left")
)
query.collect().write_csv("proteomics_map_to_strains.csv")
