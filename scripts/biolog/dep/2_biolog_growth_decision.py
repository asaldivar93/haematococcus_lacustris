import importlib.resources
import json
from pathlib import Path
from math import ceil, sqrt

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl
from plotly.subplots import make_subplots
from scipy.stats import ttest_ind
from xlsxwriter import Workbook

with importlib.resources.open_text("labutils", "plot_template.json") as file:
  style = json.load(file)

#style["layout"]["colorway"] = random.shuffle(style["layout"]["colorway"])
pio.templates["paper"] = go.layout.Template(
    data=style["data"],
    layout=style["layout"],
)
pio.templates.default = "simple_white+paper"

# By replicate
rates_df = pl.read_csv("data/1_interim/biolog/E260123/rates_df.csv")
gpb = rates_df.group_by("Plate", "Replicate", "Wavelength")
df_list = []
for name, data in gpb:
    this_data = (
        data
        .sort("Well")
        .with_columns(mu_ref=pl.col("mu (day-1)").first()) # Well A1 is the control
        .with_columns(growth=(pl.col("mu (day-1)") > pl.col("mu_ref")))
    )
    df_list.append(this_data)
rates_df = pl.concat(df_list)

# If more than 2 replicates are True for each well
thrs = 2
gpb = rates_df.group_by("Plate", "Well")
df_list = []
for name, data in gpb:
    growth_590 = sum(data.filter(pl.col("Wavelength")==590)["growth"]) > thrs
    growth_740 = sum(data.filter(pl.col("Wavelength")==740)["growth"]) > thrs
    new_dict = {
        "Plate": name[0],
        "Well": name[1],
        "growth": all([growth_590, growth_740]),
    }
    df_list.append(new_dict)
decision_sum = pl.DataFrame(df_list, orient="row")
decision_sum = decision_sum.with_columns(pl.lit("Replicates").alias("method"))
decision_sum.filter(decision_sum["growth"])

# Growth decision based on avg of replicates
rates_avg = (
    rates_df
    .group_by("Plate", "Well", "Wavelength")
    .agg(
        pl.col("mu (day-1)").mean().alias("mu_avg"),
        pl.col("mu (day-1)").std().alias("mu_std")
    )
)
gpb = rates_avg.group_by("Plate", "Wavelength")
df_list = []
for name, data in gpb:
    this_data = (
        data
        .sort("Well")
        .with_columns(
            pl.col("mu_avg").first().alias("mu_ref"),
            pl.col("mu_std").first().alias("std_ref"),
        )
        .with_columns(
            growth = pl.col("mu_avg") > (pl.col("mu_ref") + pl.col("std_ref"))
        )
    )
    df_list.append(this_data)
growth_avg = pl.concat(df_list)

decision_avg = (
    growth_avg
    .group_by("Plate", "Well")
    .agg(pl.col("growth"))
    .with_columns(
        pl.col("growth").list.all(),
        pl.lit("average").alias("method")
    )
)

# Growth decision based on t-test
p_thr = 0.05
log2_fc_thr = 0.5

gpb = rates_df.group_by("Plate", "Wavelength")
df_list = []
for name, data in gpb:
    # Get mu vals for well A1
    mu_ref = data.filter(pl.col("Well")=="A01")["mu (day-1)"].to_numpy()
    gpb_2 = data.group_by("Well")
    for well, data_2 in gpb_2:
        # Get mu values for every other well
        mu_test = data_2["mu (day-1)"].to_numpy()
        # Welch's t-test and get log2(test/ref)
        stat, p = ttest_ind(mu_ref, mu_test, equal_var=False)
        log2_fc = np.log2(np.mean(mu_test) / np.mean(mu_ref))
        # Save the data
        new_dict = {
            "Plate": name[0],
            "Well": well[0],
            "Wavelength": name[1],
            "p": p,
            "log10_p": -np.log10(p),
            "log2_fc": log2_fc,
            # growth is True if p-val and log2FC are above the thresholds
            "growth": bool((p < p_thr) and (log2_fc > log2_fc_thr))
        }
        df_list.append(new_dict)
growth_pval = pl.DataFrame(df_list, orient="row")
decision_pval = (
    growth_pval
    .group_by("Plate", "Well")
    .agg(pl.col("growth"))
    .with_columns(
        pl.col("growth").list.all(),
        pl.lit("pval").alias("method")
    )
)

# Bring all decision from the different methods together
mixed = pl.concat([decision_sum, decision_avg, decision_pval])
decision_all = (
    mixed
    .group_by("Plate", "Well")
    .agg(pl.col("growth"))
    .with_columns(
        pl.col("growth").list.all(),
        pl.lit("all").alias("method")
    )
)
decision_all.filter(decision_all["growth"]).sort("Plate", "Well")

biolog_map = pl.read_excel("data/external/biolog_map.xlsx")
q = (
    decision_all.filter(decision_all["growth"])
    .join(biolog_map, on=["Plate", "Well"], how="left")
    .sort("Plate", "Well")
    .drop("well")
)
q["growth"].sum()

out_dir = Path("data/3_results/biolog/E260123/")
out_dir.mkdir(parents=True, exist_ok=True)
file_path = out_dir / "biolog_decision.xlsx"
with Workbook(file_path) as wb:
    rates_df.write_excel(wb, worksheet="rates")
    growth_avg.write_excel(wb, worksheet="average")
    growth_pval.write_excel(wb, worksheet="pval")
    q.write_excel(wb, worksheet="decision")


out_dir = Path("data/3_results/biolog/E260123/")
corrected_od = pl.read_csv("data/1_interim/biolog/E260123/data_df.csv")
corrected_od = corrected_od.filter(pl.col("Wavelength")==590)
mixed = pl.read_excel(out_dir / "biolog_decision.xlsx", sheet_name="decision")
mixed

n_cols = ceil(sqrt(len(mixed)))
query = corrected_od.join(mixed, on=["Plate", "Well"]).with_columns(label=pl.col("Plate")+pl.lit(":")+pl.col("metabolite")).sort("metabolite")
titles = query["label"].unique().to_list()
titles.sort()
fig = make_subplots(rows=n_cols, cols=n_cols, subplot_titles=titles)
gpb = query.group_by("metabolite", maintain_order=True)
n = 1
row = 1
col = 1
for name, data in gpb:
    x=data["Time"]
    y=data["abs"]
    fig.add_trace(
        go.Scatter(
            x=x, y=y, mode="markers",
            showlegend=False
        ),
        row=row, col=col
    )
    n += 1
    if n>n_cols:
        row += 1
        n=1
    col += 1
    if col>n_cols:
        col=1
fig.update_layout({
    "height": 280 * (600 / 158.75),
    "width": 280 * (600 / 158.75),
})
fig.show()
fig.write_image(out_dir / "pos_mets.png")
fig.write_html(out_dir / "pos_mets.html")
