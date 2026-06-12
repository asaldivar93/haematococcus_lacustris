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

def plot_fit(fig, data, row, col):
    line_color = "#80b1d3"
    fill_color = "gray"
    if data["growth"][0]:
        line_color = "#33a02c"
        fill_color = "#b2df8a"

    fig.add_trace(
        go.Scatter(
            x=data_3[time],
            y=data_3["y_pred"],
            line={"color": line_color},
            mode="lines",
            showlegend=False,
        ),
        row=row, col=col,
    )

    x_fill = data_3[time].to_list() + data_3[time][::-1].to_list()
    y_upper = (data_3["y_pred"] + data_3["error"]).to_list()
    y_lower = (data_3["y_pred"] - data_3["error"])[::-1].to_list()
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
            x=data_3[time],
            y=data_3["y_log"],
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

if __name__=="__main__":
    data_file = Path("data/2_processed/E260123_biolog_light/growth_decision.xlsx")
    data_df = pl.read_excel(data_file, sheet_name="data")
    results_df = pl.read_excel(data_file, sheet_name="all")
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
    time="days"
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
            # Plot each replicate individually
            gpb_3 = data_2.group_by("replicate", maintain_order=True)
            for name_3, data_3 in gpb_3:
                fig.update_yaxes(range=[-0.5, 5], row=row, col=col)
                # Plot well A1 - Control
                control_df = (
                    data_1
                    .filter(
                        pl.col("well")=="A1",
                        pl.col("replicate")==name_3[0]
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
                fig = plot_fit(fig, data_3, row, col)
            counter += 1
            col += 1
            if col > 6:
                col = 1
                row += 1
            if counter > 31:
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

    for plate, fig_list in figures_dict.items():
        for i, fig in enumerate(fig_list):
            fig.write_image(f"data/2_processed/E260123_biolog_light/{plate}_{i}.svg", scale=3)
            fig.write_image(f"data/2_processed/E260123_biolog_light/{plate}_{i}.png", scale=3)
            fig.write_html(f"data/2_processed/E260123_biolog_light/{plate}_{i}.html")
