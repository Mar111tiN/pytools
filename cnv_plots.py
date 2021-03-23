import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# helper functions for plot_cnv
from plot_helpers import extract_pos, sort_df, get_chrom_df, make_color_chroms


def plot_genomic(
    df,
    plots,
    chroms="all",
    color_chroms=True,
    colormap="coolwarm_r",
    region="",
    figsize=(20, 4),
    ylim=(-1, 1),
    label_size=12,
):

    #### DATA MANGELING ##########
    # get cols for rearranging
    org_cols = list(df.columns)

    # sort the df
    df = sort_df(df)
    # reduce the df to the selected chromosomes
    if region:
        chrom, start, end = extract_pos(region)
        df = df.query("Chr == @chrom and @start <= Pos <= @end")
    elif chroms != "all":
        df = df.query("Chr in @chroms")

    # get the chrom_df for collapsing the
    chrom_df = get_chrom_df(df)

    df = df.merge(chrom_df.loc[:, "dif"], on="Chr")
    df["PlotPos"] = df["FullExonPos"] - df["dif"]

    # rearrange the df as return value
    new_cols = org_cols[:4] + ["PlotPos"] + org_cols[4:]
    df = df.loc[:, new_cols]

    ######## PLOTTING #######
    # plot the figure
    fig, ax = plt.subplots(figsize=figsize)

    # set the x-axis limits
    _ = ax.set_xlim(0, df["PlotPos"].max())

    # plot the graphs #######
    for plot in plots:
        if plot["plot_type"] == "line":
            plot = ax.plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
        elif plot["plot_type"] == "scatter":
            plot = ax.scatter(df["PlotPos"], df[plot["data"]], **plot["plot_args"])

    _ = ax.set_ylim(ylim)
    # add the color chroms
    _ = make_color_chroms(
        ax, chrom_df, color_chroms, ylimits=ax.get_ylim(), colormap=colormap
    )

    ######## LABELS ###################
    # set the axis labels
    _ = ax.set_xlabel("genomic coords", fontsize=1.25 * label_size)
    # quick fix for one y-label
    _ = ax.set_ylabel(
        " / ".join([plot["title"] for plot in plots]), fontsize=1.25 * label_size
    )

    # ####### CHROM LABELS #############
    add_chrom_labels(ax, chrom_df, ax.get_ylim())

    # ###### X-AXIS ####################
    # set major ticks and grid for chrom

    ax = set_ticks(ax, df, chrom_df, label_size=label_size)

    # return fig and ax for further plotting and return edited dataframe
    return fig, ax, df, chrom_df


def plot_snp(
    df,
    snp_plots=[],
    cov_plots=[],
    chroms="all",
    cov_offset=0.25,
    cov_height=0.5,
    color_chroms=True,
    colormap="coolwarm_r",
    region="",
    label_size=12,
    figsize=(20, 4),
    ylim=(-1, 1),
    fig=None,
):

    MAXLOG2RATIO = 2.5
    #### DATA MANGELING ##########
    # get cols for rearranging
    org_cols = list(df.columns)

    # sort the df
    df = sort_df(df)
    # reduce the df to the selected chromosomes
    if region:
        chrom, start, end = extract_pos(region)
        df = df.query("Chr == @chrom and @start <= Pos <= @end")
    elif chroms != "all":
        df = df.query("Chr in @chroms")

    # get the chrom_df for collapsing the
    chrom_df = get_chrom_df(df)
    df = df.merge(chrom_df.loc[:, "dif"], on="Chr")
    df["PlotPos"] = df["FullExonPos"] - df["dif"]
    # rearrange the df as return value
    new_cols = org_cols[:4] + ["PlotPos"] + org_cols[4:]
    df = df.loc[:, new_cols]

    #########################
    # ####### PLOTTING #######
    # plot the figure
    if not fig:
        fig, ax = plt.subplots(figsize=figsize)

    # set the x-axis limits
    _ = ax.set_xlim(0, df["PlotPos"].max())

    # PLOT COV Data
    if len(cov_plots):
        scale_factor = cov_height / (MAXLOG2RATIO + 1)
        offset = 1 + scale_factor + cov_offset

        ylim = (ylim[0], ylim[1] + cov_offset + cov_height)

        for plot in cov_plots:
            # normalize the coverage data:
            # 2.5 is the approx max log2ratio (LOH to 8N)

            df[plot["data"]] = df[plot["data"]] * scale_factor + offset
            if plot["plot_type"] == "line":
                plot = ax.plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
            elif plot["plot_type"] == "scatter":
                # highjack plot_args
                pa = plot["plot_args"]
                if "c" in pa:
                    pa["c"] = df[pa["c"]]
                if "s" in pa:
                    if isinstance(pa["s"], str):
                        pa["s"] = df[pa["s"]] * 20 + 1
                plot = ax.scatter(df["PlotPos"], df[plot["data"]], **pa)

    ######## plot the SNP graphs #######
    for plot in snp_plots:
        if plot["plot_type"] == "line":
            plot = ax.plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
        elif plot["plot_type"] == "scatter":
            # highjack plot_args with
            pa = plot["plot_args"]
            if "c" in pa:
                pa["c"] = df[pa["c"]]
            if "s" in pa:
                if isinstance(pa["s"], str):
                    pa["s"] = df[pa["s"]] * 20 + 1
            plot = ax.scatter(df["PlotPos"], df[plot["data"]], **pa)

    _ = ax.set_ylim(ylim)
    # add the color chroms
    _ = make_color_chroms(
        ax, chrom_df, color_chroms, ylimits=ax.get_ylim(), colormap=colormap
    )

    # ####### LABELS ###################
    # set the axis labels
    _ = ax.set_xlabel("genomic coords", fontsize=1.25 * label_size)
    # quick fix for one y-label
    _ = ax.set_ylabel(
        " / ".join([plot["title"] for plot in snp_plots]), fontsize=1.25 * label_size
    )

    # ####### CHROM LABELS #############
    add_chrom_labels(ax, chrom_df, ax.get_ylim())

    # ###### X-AXIS ####################
    # set major ticks and grid for chrom

    ax = set_ticks(ax, df, chrom_df, label_size=label_size)

    # return fig and ax for further plotting and return edited dataframe
    return fig, ax, df, chrom_df


def plot_snp2(
    df,
    snp_plots=[],
    cov_plots=[],
    blocks=[],
    chroms="all",
    cov_offset=0.25,
    cov_height=0.5,
    color_chroms=True,
    colormap="coolwarm_r",
    region="",
    label_size=12,
    figsize=(20, 4),
    ylim=(-1, 1),
):

    MAXLOG2RATIO = 2.5
    #### DATA MANGELING ##########
    # get cols for rearranging
    org_cols = list(df.columns)

    # sort the df
    df = sort_df(df)
    # reduce the df to the selected chromosomes
    if region:
        chrom, start, end = extract_pos(region)
        df = df.query("Chr == @chrom and @start <= Pos <= @end")
    elif chroms != "all":
        df = df.query("Chr in @chroms")

    # get the chrom_df for collapsing the
    chrom_df = get_chrom_df(df)
    df = df.merge(chrom_df.loc[:, "dif"], on="Chr")
    df["PlotPos"] = df["FullExonPos"] - df["dif"]
    # rearrange the df as return value
    new_cols = org_cols[:4] + ["PlotPos"] + org_cols[4:]
    df = df.loc[:, new_cols]

    #########################
    ######## PLOTTING #######
    # plot the figure
    fig, ax = plt.subplots(figsize=figsize)

    # set the x-axis limits
    _ = ax.set_xlim(0, df["PlotPos"].max())

    # PLOT COV Data
    if len(cov_plots):
        scale_factor = cov_height / (MAXLOG2RATIO + 1)
        offset = 1 + scale_factor + cov_offset

        ylim = (ylim[0], ylim[1] + cov_offset + cov_height)

        for plot in cov_plots:
            # normalize the coverage data:
            # 2.5 is the approx max log2ratio (LOH to 8N)

            df[plot["data"]] = df[plot["data"]] * scale_factor + offset

            minus, zero, one = [c * scale_factor + offset for c in [-1, 0, 1]]
            if plot["plot_type"] == "line":
                plot = ax.plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
            elif plot["plot_type"] == "scatter":
                # highjack plot_args
                pa = plot["plot_args"]
                if "c" in pa:
                    pa["c"] = df[pa["c"]]
                if "s" in pa:
                    if isinstance(pa["s"], str):
                        pa["s"] = df[pa["s"]] * 20 + 1
                plot = ax.scatter(df["PlotPos"], df[plot["data"]], **pa)

    ######## plot the SNP graphs #######
    for plot in snp_plots:
        if plot["plot_type"] == "line":
            plot = ax.plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
        elif plot["plot_type"] == "scatter":
            # highjack plot_args with
            pa = plot["plot_args"]
            if "c" in pa:
                pa["c"] = df[pa["c"]]
            if "s" in pa:
                if isinstance(pa["s"], str):
                    pa["s"] = df[pa["s"]] * 20 + 1
            plot = ax.scatter(df["PlotPos"], df[plot["data"]], **pa)

    _ = ax.set_ylim(ylim)
    # add the color chroms
    _ = make_color_chroms(
        ax, chrom_df, color_chroms, ylimits=ax.get_ylim(), colormap=colormap
    )

    ######## LABELS ###################
    # set the axis labels
    _ = ax.set_xlabel("genomic coords", fontsize=1.25 * label_size)
    # quick fix for one y-label
    _ = ax.set_ylabel(
        " / ".join([plot["title"] for plot in snp_plots]), fontsize=1.25 * label_size
    )

    ######## BLOCKS ##################
    if len(blocks):
        for col in blocks:
            if col.startswith("snp"):
                df.loc[df[col] > 0, col] = 0.5
                df.loc[df[col] == 0, col] = -10
                plot = ax.scatter(df["PlotPos"], df[col], s=15, color="green")
        if col.startswith("cov"):
            df.loc[df[col] > 0, col] = offset
            df.loc[df[col] == 0, col] = -10
            plot = ax.scatter(df["PlotPos"], df[col], s=10, color="blue")

    ######## CHROM LABELS #############
    add_chrom_labels(ax, chrom_df, ax.get_ylim())

    ####### X-AXIS ####################
    # set major ticks and grid for chrom

    ax = set_ticks(ax, df, chrom_df, label_size=label_size)
    # set helper lines
    _ = ax.axhline(y=1, c="k", lw=2, ls="-")
    _ = ax.axhline(y=0.5, c="k", lw=1.5, alpha=0.5, ls="--")

    _ = ax.axhline(y=minus, c="k", lw=1.5, alpha=0.5, ls="--")
    _ = ax.axhline(y=zero, c="k", lw=1.5, alpha=0.5, ls="-")
    _ = ax.axhline(y=one, c="k", lw=1.5, alpha=0.5, ls="--")
    # return fig and ax for further plotting and return edited dataframe
    return fig, ax, df, chrom_df


def plot_2d(df, xcol, ycol, df2=pd.DataFrame(), figsize=(5, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    _ = ax.scatter(df[xcol], df[ycol], s=0.1)
    if len(df2.index):
        _ = ax.scatter(df2[xcol], df2[ycol], s=1, alpha=0.5, color="red")
    _ = ax.set_xlabel(xcol, fontsize=10)
    _ = ax.set_ylabel(ycol, fontsize=10)

    def get_lims(col):
        if "log" in col:
            return (-1.5, 3)
        if "abs" in col:
            return (0, 1)
        if col == "deltaVAFvar":
            return (0, 0.2)
        if col == "deltaVAFstd":
            return (0, 1)
        if col == "VAF":
            return (0, 1)
        else:
            return (-1, 1)

    _ = ax.set_xlim(get_lims(xcol))
    _ = ax.set_ylim(get_lims(ycol))


def plot_3d(df, xcol, ycol, zcol, df2=pd.DataFrame(), figsize=(10, 10)):
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    _ = ax.scatter3D(df[xcol], df[ycol], df[zcol], color="green", alpha=0.2, s=0.1)
    if len(df2.index):
        _ = ax.scatter3D(df2[xcol], df2[ycol], df2[zcol], s=5, color="red")
    # labels
    _ = ax.set_xlabel(xcol, fontsize=10)
    _ = ax.set_ylabel(ycol, fontsize=10)
    _ = ax.set_zlabel(zcol, fontsize=10)

    def get_lims(col):
        if "log" in col:
            return (-1.5, 3)
        if "abs" in col:
            return (0, 1)
        if col == "deltaVAFvar":
            return (0, 0.2)
        if col == "deltaVAFstd":
            return (0, 1)
        else:
            return (-1, 1)

    _ = ax.set_xlim(get_lims(xcol))
    _ = ax.set_ylim(get_lims(ycol))
    _ = ax.set_zlim(get_lims(zcol))
    return fig, ax
