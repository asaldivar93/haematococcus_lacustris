"""Annotation based on best hit to swissprot.

TO DO:
    - replace to orthologues matching to eggnog?
    - do orthogroups in swissprot and replace to orthomatch to that?
"""

import subprocess
from pathlib import Path

import polars as pl

cdd_map_file = Path.home() / "databases/cdd/cddid.tbl"
CDD_DB = pl.scan_csv(
    cdd_map_file,
    separator="\t",
    has_header=False,
    new_columns=["number", "target", "name", "description", "number2"],
).select("target", "name", "description")

SWP_DB = pl.scan_parquet(
    Path.home() / "databases/swissprot/uniprot_map/*.parquet",
)


def build_search_cmd(
    query_db,
    target_db,
    results_db,
    sensitivity=7.5,
    threads=12,
):
    cmd = [
        "mmseqs",
        "search",
        query_db,
        target_db,
        results_db,
        "tmp",
        "-s",
        str(sensitivity),
        "--threads",
        str(threads),
    ]
    return " ".join(cmd)


def build_convertails_cmd(query_db, target_db, results_db, out_file, out_fmt):
    cmd = [
        "mmseqs",
        "convertalis",
        query_db,
        target_db,
        results_db,
        out_file,
        "--format-output",
        out_fmt,
    ]
    return " ".join(cmd)


def run_command_realtime(cmd):
    """Run a command and prints its output in real time."""
    # Use Popen to control the process flow
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Redirect stderr to stdout to capture all output
        text=True,  # Decode output as text (Python 3.6+)
        bufsize=1,  # Line bufferring
        shell=True,
    )

    # Read output line by line as it becomes available
    stdout = ""
    for line in process.stdout:
        stdout += line
        print(
            line,
            end="",
            flush=True,
        )  # Print immediately with no extra newlines

    # Wait for the process to finish and get the exit code
    process.wait()
    return process, stdout


if __name__ == "__main__":
    # genomes to annotate
    genomes = pl.read_csv("data/genomes.list", separator="\t")
    directories = {
        "interim": Path("data/1_interim/annotation"),
        "results": Path("data/2_results/annotation"),
    }

    for path in directories.values():
        path.mkdir(exist_ok=True, parents=True)

    # Create a mmseqs database if it does not exist already
    db_file = directories["interim"] / "query.db"
    if not db_file.is_file():
        query_files = " ".join(genomes["path"].to_list())
        cmd = ["mmseqs", "createdb", query_files, str(db_file)]
        cmd = " ".join(cmd)
        result, stdout = run_command_realtime(cmd)

    # Annotate on swissprot database
    # Get arguments
    format_out = "query,target,pident,evalue,bits,theader,qset,qlen,tlen,qcov,tcov,qstart,qend,tstart,tend"
    query_db = str(db_file)
    target_db = "$HOME'/databases/swissprot/swissprot'"
    results_db = str(directories["interim"] / "swissprot_all_hits.db")
    results_file = str(directories["interim"] / "swissprot_all_hits.tsv")
    # Run annotation
    cmd = build_search_cmd(
        query_db,
        target_db,
        results_db,
        sensitivity=7.5,
        threads=n_cores,
    )
    result, stdout = run_command_realtime(cmd)
    if result.returncode != 0:
        raise ValueError(f"Error in {cmd}")
    # Convert to tsv
    cmd = build_convertails_cmd(
        query_db,
        target_db,
        results_db,
        results_file,
        format_out,
    )
    result, stdout = run_command_realtime(cmd)
    if result.returncode != 0:
        raise ValueError(f"Error in {cmd}")

    # Annotate on cdd database
    # Get arguments
    target_db = "$HOME'/databases/cdd/cdd.db'"
    results_db = str(directories["interim"] / "cdd_all_hits.db")
    results_file = str(directories["interim"] / "cdd_all_hits.tsv")
    # Run Annotation
    cmd = build_search_cmd(
        query_db,
        target_db,
        results_db,
        sensitivity=7.5,
        threads=28,
    )
    result, stdout = run_command_realtime(cmd)
    if result.returncode != 0:
        raise ValueError
    # Convert to tsv
    cmd = build_convertails_cmd(
        query_db,
        target_db,
        results_db,
        results_file,
        format_out,
    )
    result, stdout = run_command_realtime(cmd)
    if result.returncode != 0:
        raise ValueError

    # Map of original file names to strain id
    strain_map = {
        row["path"].split("/")[-1]: row["strain"]
        for row in genomes.iter_rows(named=True)
    }

    # Add CDD annotation
    results_file = directories["interim"] / "cdd_all_hits.tsv"
    cdd_hits_lazy = pl.scan_csv(
        results_file,
        separator="\t",
        has_header=False,
        new_columns=format_out.split(","),
    )
    query = (
        cdd_hits_lazy
        # Replace file name with strain id
        .with_columns(pl.col("qset").replace(strain_map))
        # Keep only top 5 hits
        .sort("query", "bits", descending=True)
        .group_by("query")
        .head(5)
        # Add CDD info
        .join(CDD_DB, on="target", how="left")
    )
    # Save one file per strain
    for s in strain_map.values():
        out_path = directories["results"] / f"cdd_{s}_hits.csv"
        query.filter(pl.col("qset") == s).sink_csv(out_path)

    # Add swissprot annotation
    results_file = str(directories["interim"] / "swissprot_all_hits.tsv")
    swissprot_hits_lazy = pl.scan_csv(
        results_file,
        separator="\t",
        has_header=False,
        new_columns=format_out.split(","),
    )

    query = (
        swissprot_hits_lazy
        # Replace file name with strain id
        .with_columns(pl.col("qset").replace(strain_map))
        # Get the top 5 hits
        .sort("query", "bits", descending=True)
        .group_by("query")
        .head(5)
        # Add swissprot info
        .join(SWP_DB, on="target", how="left")
        .with_columns(
            pl.col("ec-code").list.join(";"),
            pl.col("location").list.join(";"),
        )
    )
    # Save one file per strain
    for s in strain_map.values():
        out_path = directories["results"] / f"swp_{s}_hits.csv"
        query.filter(pl.col("qset") == s).sink_csv(out_path)
