import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
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


##########################################################################
############ plot GISTIC like CNV plots ##################################

def add_chrom_blocks(ax, chrom_df, colormap="coolwarm_r", ylimits=(-10, 10), yoffset=5, label_size=1):
    # set the cmap from provided argument
    cmap = plt.cm.get_cmap("coolwarm_r", 23)
    rects = []
    ylimits=ax.get_ylim()
    ymin = ylimits[0] * 1.1
    height = (ylimits[1] - ymin) * 1.1
    for _, row in chrom_df.iterrows():
        rect = Rectangle(
            (row["FullStart"], ymin), width=row["FullEnd"] - row["FullStart"], height=height
            )
        rects.append(rect)
    rect_kwargs = dict(alpha=1, fc="none", ec="black", lw=1, ls="dotted")
    rect_collection = PatchCollection(rects, **rect_kwargs)
    # rect_collection.set_array(
     #    chrom_df['Chr'].str.replace("chr", "").str.replace("X", "23").astype(int)
    # )
    _ = ax.add_collection(rect_collection)
    
    
    # add the labels
    # get the min_chrom_fraction from minimum chrom_size / whole stretch
    min_chrom_frac = (chrom_df["FullEnd"] - chrom_df["FullStart"]).min() / chrom_df["FullEnd"].max()
    chrom_size = min(20, max(15, label_size * 200 * min_chrom_frac))
    
    label_style = dict(size=chrom_size, color="#2f3832")
    # set the height and ymin beyond the ylimits so borders are not seen
    ypos = ylimits[0] - yoffset
    for _, row in chrom_df.iterrows():
        mean = (row['FullEnd'] + row['FullStart']) / 2
        chrom = row['Chr'].replace("chr", "")
        ax.text(mean, ypos, chrom, ha="center", **label_style)   
    
    return ax


def add_cytoband_blocks(ax, band_df, yoffset=5, label_size=1, band_colors={"p":'lightgray', "q":'darkgray'}):
    # set the cmap from provided argument
    # cmap = plt.cm.get_cmap("coolwarm_r", 23)
    
    ylimits=ax.get_ylim()
    ymin = ylimits[0] * 1.1
    height = (ylimits[1] - ymin) * 1.1
    ypos = ylimits[0] + yoffset
    
    label_style = dict(size=10 * label_size, color="#2f3832")
    for band in ['p', 'q']:
        rects = []
        for _, row in band_df.query('BAND == @band').iterrows():
            rect = Rectangle(
                (row["FullStart"], ymin), width=row["FullEnd"] - row["FullStart"], height=height
                )
            rects.append(rect)
            # add the label
            mean = (row['FullEnd'] + row['FullStart']) / 2
            ax.text(mean, ypos, band, ha="center", **label_style)
        
        # each rect collection has separate rect_kwargs
        rect_kwargs = dict(alpha=.4, fc=band_colors[band], ec="none", lw=.3, ls="-")
        rect_collection = PatchCollection(rects, **rect_kwargs)
        _ = ax.add_collection(rect_collection)
    
    return ax


def get_tick_range(df, col):
    _max = df[col].max()
    # positive values
    if _max > 0:
        step = 10 if _max >= 20 else 5

        return [s + step for s in range(0,step * math.ceil(_max / step),step)]
    else:
        _min = df[col].min()
        step = 10 if _min <= -20 else 5
        return [s for s in range(step * math.floor(_min / step),0,step)]

    
def set_y_ticks(ax, df, label_size=12, margin=2):
    """
    for a given tick number, set nicely spread ticks
    """

    # get the ticks from the gains and losses range
    y_ticks = get_tick_range(df, "losses") + get_tick_range(df, "gains")
    
    # set the ylim from y_ticks range
    _ = ax.set_ylim(y_ticks[0] - margin, y_ticks[-1] + margin)

    ax.yaxis.set_major_locator(plt.FixedLocator(y_ticks))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(np.abs(y_ticks)))
    
    ax.yaxis.grid(which="major", linestyle="dotted", linewidth=1)

    ax.yaxis.set_tick_params(which="major", length=20, labelsize=label_size)
    # set the tick labels
    for tick in ax.yaxis.get_majorticklabels():
        tick.set_horizontalalignment("left")
    return ax

        
def make_CNV_plot(df, chrom_df, band_df,
                  figsize=(10,4),
                  y_margin=2,
                  label_size=10,
                  cnv_fill_colors={},
                  chrom_color_scheme="coolwarm_r",
                  chrom_label_yoffset=5,
                  chrom_label_size=1,
                  band_label_yoffset=5,
                  band_label_size=1,
                  band_colors={"p":'lightgray', "q":'darkgray'},
                  y_label_size=10
                 ):
    '''
    
    '''
    
    ######### data ############
    # get the FullPos from chrom_df into df
    df['FullPos'] = df['Start'] + df.merge(chrom_df, on="Chr")['FullStart'] - 1
    
    # get the FullPos from chrom_df into band_df
    band_df = band_df.merge(chrom_df, on="Chr")
    band_df['FullEnd'] = band_df['End'] + band_df['FullStart'] - 1
    band_df['FullStart'] = band_df['Start'] + band_df['FullStart']
    band_df = band_df.loc[:, ['Chr', 'BAND', 'FullStart', 'FullEnd']]
    
    # create the figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # remove seaborn grid
    ax.grid(False)
    # set the x_lims
    _ = ax.set_xlim(0, chrom_df['FullEnd'].max())
    
    _ = set_y_ticks(ax, df, label_size=y_label_size, margin=y_margin)
    
    # plot the graphs
    plot = ax.stackplot(df['FullPos'], df['gains'], alpha=0.9, color=cnv_fill_colors['gains'])
    plot = ax.stackplot(df['FullPos'], df['losses'], alpha=0.9, color=cnv_fill_colors['losses']) 
        
    _ = add_chrom_blocks(ax, chrom_df, 
                         colormap=chrom_color_scheme, 
                         ylimits=ax.get_ylim(), 
                         yoffset=chrom_label_yoffset, 
                         label_size=chrom_label_size
                        )
    
    _ = add_cytoband_blocks(ax, band_df,
                            yoffset=band_label_yoffset,
                            label_size=band_label_size,
                            band_colors=band_colors
                           )
    

    
    # remove the x-plot ticks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)
    
    
    return fig, ax
