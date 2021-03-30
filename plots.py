import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# helper functions for plot_cnv
from plot_helpers import *


######### 2D-plots for SNP/VAF ##################################
def make_SNP_plot(file1, file2, plot_file="", plot_quality=90):
    """
    plot the difference in SNPs between two samples and calculates offRate
    """

    # load the tumor file
    t_vaf = pd.read_csv(file1, sep="\t").loc[:, ["FullExonPos", "VAF"]]
    n_vaf = pd.read_csv(file2, sep="\t").loc[:, ["FullExonPos", "VAF"]]
    # merge for corresponding SNP pos
    t_n = t_vaf.merge(n_vaf, on="FullExonPos").drop("FullExonPos", axis=1)

    sample_t = os.path.basename(file1).replace(".snp", "").replace("_", "")
    sample_n = os.path.basename(file2).replace(".snp", "").replace("_", "")

    fig, ax = plt.subplots(figsize=(10, 10))
    _ = ax.scatter(t_n["VAF_x"], t_n["VAF_y"], s=0.2, alpha=0.2)
    _ = ax.set_xlabel(sample_t, fontsize=20)
    _ = ax.set_ylabel(sample_n, fontsize=20)
    # calculate offRate
    df0 = t_n[(t_n > 0.1).any(axis=1)]
    n = len(df0.index)
    df1 = df0[np.abs(df0["VAF_x"] - df0["VAF_y"]) > 0.25]
    m = len(df1.index)
    off_ratio = m / n * 100
    _ = ax.set_title(
        f"{sample_t} vs {sample_n} - offRate {round(off_ratio, 1)}", fontsize=30
    )

    # plot the graph
    if plot_file:
        fig.savefig(plot_file, quality=plot_quality)
    return fig


def plot_VAFs(df, sample="", plot_file="", plot_quality=90):
    """
    takes a dataframe and
    """
    if sample:
        df = df.query("sample == @sample")
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 10))
    _ = ax.scatter(df["TVAF"], df["NVAF"], s=0.5, alpha=0.5)
    _ = ax.set_xlabel("TVAF", fontsize=20)
    _ = ax.set_ylabel("NVAF", fontsize=20)
    _ = ax.set_title(sample, fontsize=30)
    if plot_file:
        fig.savefig(plot_file, quality=plot_quality)
    return fig













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




