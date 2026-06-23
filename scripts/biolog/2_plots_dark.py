"""Plots of biolog data."""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl
from plotly.subplots import make_subplots

with Path("plot_template.json").open("r") as file:
    style = json.load(file)

pio.templates["paper"] = go.layout.Template(
    data=style["data"],
    layout=style["layout"],
)
pio.templates.default = "simple_white+paper"
pio.renderers.default = "browser"

def plot_fit(fig, data, row, col, time="days"):
    line_color = "#80b1d3"
    fill_color = "gray"
    if data["growth"][0]:
        line_color = "#33a02c"
        fill_color = "#b2df8a"

    fig.add_trace(
        go.Scatter(
            x=data[time],
            y=data["y_pred"],
            line={"color": line_color},
            mode="lines",
            showlegend=False,
        ),
        row=row, col=col,
    )

    x_fill = data[time].to_list() + data[time][::-1].to_list()
    y_upper = (data["y_pred"] + data["error"]).to_list()
    y_lower = (data["y_pred"] - data["error"])[::-1].to_list()
    y_fill = y_upper + y_lower

    fig.add_trace(
        go.Scatter(
            x=x_fill, y=y_fill,
            fill="toself", opacity=0.5,
            line={"color": fill_color, "dash": "dot"},
            mode="lines",
            showlegend=False,
        ),
        row=row, col=col
    )

    fig.add_trace(
        go.Scatter(
            x=data[time],
            y=data["y_log"],
            line={"color": "#fb8072"},
            marker={"size": 4},
            mode="markers",
            showlegend=False,
        ),
        row=row, col=col,
    )

    well = data["well"][0]
    met_name = data["metabolite"][0]
    name = f"{well}:{met_name}" if len(met_name)<20 else well
    fig.add_annotation(
        text=name,
        x=5, y=4,
        showarrow=False,
        font={"size": 8},
        row=row, col=col
    )
    return fig

def plot_biolog_data(data_df, results_df, time="days"):
    # Add growth decision
    query = (
        data_df
        .join(
            results_df.select("plate", "well", "growth"),
            on=["plate", "well"],
        )
        .sort("plate")
    )

    # Plot growth curves
    gpb_1 = query.group_by("plate", maintain_order=True)
    figures_dict = {}
    col=1
    row=1
    counter=0
    fig = make_subplots(rows=6, cols=6)
    for name_1, data_1 in gpb_1:
        # Split wells in 3 plots of 32 wells each for each plate
        figures_dict[name_1[0]] = []
        gpb_2 = data_1.sort("well").group_by("well", maintain_order=True)
        for name_2, data_2 in gpb_2:
            fig.update_yaxes(range=[-0.5, 5], row=row, col=col)
            # Plot well A1 - Control
            control_df = (
                data_1
                .filter(
                    pl.col("well")=="A01",
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=control_df[time],
                    y=control_df["y_pred"],
                    line={"color": "#cab2d6"},
                    mode="lines",
                    showlegend=False,
                ),
                row=row, col=col,
            )
            # Plot each Well
            fig = plot_fit(fig, data_2, row, col)
            counter += 1
            col += 1
            if col > 5:
                col = 1
                row += 1
            if counter > 23:
                counter = 0
                row = 1
                col = 1
                fig.update_layout({
                    "height": 175 * (600 / 158.75),
                    "width": 175 * (600 / 158.75),
                    "font": {"size": 10},
                })
                figures_dict[name_1[0]].append(fig)
                fig = make_subplots(rows=6, cols=6)

    return figures_dict

def get_size_from_mm(height_mm, width_mm):
    mm_coversion = (600 / 158.75)
    return height_mm * mm_coversion, width_mm * mm_coversion


def plot_volcano(p_val_df):
    fig = px.scatter(
        p_val_df.filter(pl.col("growth")),
        x="log2_fc",
        y="log10_p",
        color="plate",
        hover_name="metabolite",
    )
    fig.add_trace(
        go.Scatter(
            x=p_val_df.filter(~pl.col("growth"))["log2_fc"],
            y=p_val_df.filter(~pl.col("growth"))["log10_p"],
            mode="markers",
            marker={"color": "gray"},
            opacity=0.5,
            showlegend=False,
        )
    )
    fig.add_vline(
        x=p_val_df.filter(pl.col("growth"))["log2_fc"].min(),
        line_width=2, line_dash="dash", line_color="black",
        opacity=1
    )
    fig.add_hline(
        y=p_val_df.filter(pl.col("growth"))["log10_p"].min(),
        line_width=2, line_dash="dash", line_color="black",
        opacity=1
    )

    return fig


if __name__=="__main__":
    input_dir = "data/2_processed/E260123_biolog_dark"
    output_dir = "data/2_processed/E260123_biolog_dark"
    input_file = "growth_decision.xlsx"
    directory = {
        "input": Path(input_dir),
        "output": Path(output_dir),
    }
    data_file = directory["input"] / input_file

    # Load the data
    data_df = pl.read_excel(data_file, sheet_name="data")
    results_df = pl.read_excel(data_file, sheet_name="all")
    p_val_df = pl.read_excel(data_file, sheet_name="pval")
    fits_df = pl.read_excel(data_file, sheet_name="growth")

    # Plot growth curves and fits
    figures_dict = plot_biolog_data(data_df, results_df)
    for plate, fig_list in figures_dict.items():
        for i, fig in enumerate(fig_list):
            fig.write_image(directory["output"] / f"{plate}_{i}.svg", scale=3)
            fig.write_image(directory["output"] / f"{plate}_{i}.png", scale=3)
            fig.write_html(directory["output"] / f"{plate}_{i}.html")

    # Plot volcano plots
    height, width = get_size_from_mm(80, 85)
    fig_pval = plot_volcano(p_val_df)
    fig_pval.update_layout({
        "xaxis": {"range": (-8, 8), "dtick": 2.5},
        "showlegend": False,
        "height": height, "width": width,
    })
    fig_pval.write_image(directory["output"] / "volcano_plot.svg", scale=3)
    fig_pval.write_image(directory["output"] / "volcano_plot.png", scale=3)
    fig_pval.write_html(directory["output"] / "volcano_plot.html")

    # Box Plot growth rates and R2
    fig_mu = px.box(
        fits_df.drop_nulls().filter(pl.col("mu_err")<1e3, pl.col("mu")<1e3),
        y="mu", x="plate", color="plate"
    )
    fig_mu.update_layout({
        "xaxis": {"title": {"text": ""}},
        "showlegend": False,
        "height": height, "width": width,
    })
    fig_mu.write_image(directory["output"] / "box_growth_rates.svg", scale=3)
    fig_mu.write_image(directory["output"] / "box_growth_rates.png", scale=3)
    fig_mu.write_html(directory["output"] / "box_growth_rates.html")

    fig_r2 = px.box(
        fits_df.drop_nulls().filter(pl.col("mu_err")<1e3),
        y="R2", x="plate", color="plate"
    )
    fig_r2.update_layout({
        "xaxis": {"title": {"text": ""}},
        "showlegend": False,
        "height": height, "width": width,
    })
    fig_r2.write_image(directory["output"] / "box_r2.svg", scale=3)
    fig_r2.write_image(directory["output"] / "box_r2.png", scale=3)
    fig_r2.write_html(directory["output"] / "box_r2.html")
