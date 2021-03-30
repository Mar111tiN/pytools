import os
import pandas as pd

#################### GLOBALS ##############
coords_pattern = r"(?P<Chr>chr[0-9XY]+):(?P<Start>[0-9]+)-(?P<End>[0-9]+)"
chrom_list = [f"chr{i}" for i in range(23)] + ["chrX", "chrY"]

def show_columns(df):
    for i, col in enumerate(df.columns):
        print(i, col)


def cat_chrom(df):
    '''
    turns Chr columns into categorical
    '''
    df.loc[:, "Chr"] = pd.Categorical(df["Chr"], chrom_list)
    return df


def ken2mut(mut_file, sheet="result", hg38_bed="", selected_cols=[], trans_dict={}):
    """
    converts a kenichi merge file (excel) into a filter file with our formats
    if no hg38_bed file is given, ken2mut will output a hg19.bed file at the path of the input file
    """

    default_dict = {
        "id": "sample",
        "Func.refGene": "Func",
        "Gene.refGene": "Gene",
        "driver": "isDriver",
        "Merge_Func": "ExonicFunc",
        "AAChange.refGene": "AAChange",
        "genomicSuperDups": "SuperDups",
        "depth_tumor": "Tdepth",
        "variantNum_tumor": "TR2",
        "depth_normal": "Ndepth",
        "variantNum_normal": "NR2",
        "misRate_tumor": "TVAF",
        "misRate_normal": "NVAF",
        "P-value(fisher_realignment)": "FisherScore",
        "P-value(EBCall)": "EBscore",
    }

    # this is the transdict, if bases_tumor is there
    short_trans_dict = {
        "id": "sample",
        "Func.refGene": "Func",
        "Gene.refGene": "Gene",
        "driver": "isDriver",
        "Merge_Func": "ExonicFunc",
        "AAChange.refGene": "AAChange",
        "genomicSuperDups": "SuperDups",
        "P-value(fisher_realignment)": "FisherScore",
        "P-value(EBCall)": "EBscore",
    }

    # read file
    mut_df = pd.read_excel(mut_file, sheet_name=sheet, engine="openpyxl").query(
        "Chr == Chr"
    )
    # convert positions
    for c in ["Start", "End"]:
        mut_df[c] = mut_df[c].astype(int)

    # convert to bed_file, if no hg38_bed is given
    mut_df["hg19"] = (
        "chr"
        + mut_df["Chr"].astype(str)
        + ":"
        + mut_df["Start"].astype(str)
        + "-"
        + mut_df["End"].astype(str)
    )
    if hg38_bed == "":

        hg19_bed_file = os.path.splitext(mut_file)[0] + ".hg19.bed"
        mut_df.loc[:, "hg19"].to_csv(hg19_bed_file, index=False)
        print(f"hg19 bed file written to {hg19_bed_file}")
        return hg19_bed_file

    # hg38.bed is given
    # read hg38 bed

    hg38_df = pd.read_csv(hg38_bed, header=None, names=["hg38"])
    hg38_df

    # merge hg38 into mut_df and re-extract the coords
    mut_df = mut_df.merge(hg38_df, left_index=True, right_index=True)

    mut_df.loc[:, ["Chr", "Start", "End"]] = mut_df["hg38"].str.extract(coords_pattern)

    # rename cols
    if not trans_dict:
        if "bases_tumor" in mut_df.columns:
            trans_dict = short_trans_dict
            # extract the TR values from bases_tumor/normal
            TR_pattern = r"([0-9]+)[.,]([0-9]+)[.,]([0-9]+)[.,]([0-9]+)"

            tumor_cols = ["TR+", "TR2+", "TR-", "TR2-"]
            mut_df.loc[:, tumor_cols] = (
                mut_df["bases_tumor"]
                .str.extract(TR_pattern)
                .rename({i: col for i, col in enumerate(tumor_cols)}, axis=1)
            )

            normal_cols = ["NR+", "NR2+", "NR-", "NR2-"]
            mut_df.loc[:, normal_cols] = (
                mut_df["bases_normal"]
                .str.extract(TR_pattern)
                .rename({i: col for i, col in enumerate(normal_cols)}, axis=1)
            )
            for col in tumor_cols + normal_cols:
                mut_df[col] = mut_df[col].astype(int)
            mut_df.loc[:, ["Tdepth"]] = mut_df["TR+"] + mut_df["TR-"]
            mut_df.loc[:, ["TR2"]] = mut_df["TR2+"] + mut_df["TR2-"]
            mut_df.loc[:, ["TR1+"]] = mut_df["TR+"] - mut_df["TR2+"]
            mut_df.loc[:, ["TR1"]] = mut_df["TR1+"] + mut_df["TR-"] - mut_df["TR2-"]
            mut_df.loc[:, ["TVAF"]] = mut_df["TR2"] / mut_df["Tdepth"]

            mut_df.loc[:, ["Ndepth"]] = mut_df["NR+"] + mut_df["NR-"]
            mut_df.loc[:, ["NR2"]] = mut_df["NR2+"] + mut_df["NR2-"]
            mut_df.loc[:, ["NR1+"]] = mut_df["NR+"] - mut_df["NR2+"]
            mut_df.loc[:, ["NR1"]] = mut_df["NR1+"] + mut_df["NR-"] - mut_df["NR2-"]
            mut_df.loc[:, ["NVAF"]] = mut_df["NR2"] / mut_df["Ndepth"]

        else:
            trans_dict = default_dict

    mut_df = mut_df.rename(trans_dict, axis=1)

    # reduce to wanted columns and reorder as such
    if selected_cols:
        mut_df = mut_df.loc[:, selected_cols]

    # sort
    mut_df["Chr"] = pd.Categorical(
        mut_df["Chr"], [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    )
    mut_df = mut_df.sort_values(["sample", "Chr", "Start"])

    return mut_df


def make_cohort(
    folder,
    filter_type="filter1",
    output="",
    gzip=False,
    selected_columns=[],
    trans_dict={},
):
    """
    aggregates all filter files of a specified type
    option: if output is given as a path, it writes out to a file
    option: if trans_dict / selected_colums is given, cols are FIRST renamed and THEN selected/ordered
    reduces to selected columns
    """

    # check folder
    if not os.path.isdir(folder):
        print(f"Folder {folder} does not exist")
        return

    # get the filter files
    if filter_type:
        files = [file for file in list(os.walk(folder))[0][2] if filter_type in file]
        output = output.replace(".csv", f".{filter_type}.csv")
    else:
        files = [file for file in list(os.walk(folder))[0][2] if "-B.csv" in file]
    # get the samples from the file list
    replace_pattern = f"-B.{filter_type}.csv" if filter_type else "-B.csv"
    # print(replace_pattern)
    samples = [f.replace(replace_pattern, "").replace("_", "") for f in files]

    # load all the filter_dfs and include sample columns
    filter_dfs = []
    for file, sample in zip(files, samples):
        print(file)
        filter_df = pd.read_csv(os.path.join(folder, file), sep="\t")
        filter_df["sample"] = sample
        filter_dfs.append(filter_df)

    # put sample first
    filter_df = pd.concat(filter_dfs).loc[:, ["sample"] + list(filter_df.columns)[:-1]]
    # print(list(filter_df.columns)) --> copy for rename
    # reduce to important columns
    # use set intersection with existing columns

    df = filter_df.rename(trans_dict, axis=1)

    if selected_columns:
        selected_columns = [col for col in selected_columns if col in df.columns]
    else:
        selected_columns = list(df.columns)
    df = df.loc[:, selected_columns]
    # make Chr categorical and sort
    df["Chr"] = pd.Categorical(
        df["Chr"], [f"chr{chrom}" for chrom in list(range(1, 23)) + ["X", "Y"]]
    )
    df = df.sort_values(["sample", "Chr", "Start", "End"])
    if output:
        if gzip:
            output = f"{output}.gz"
            df.to_csv(output, sep="\t", index=False, compression="gzip")
        else:
            df.to_csv(output, sep="\t", index=False)
        print(f"Written mutation file to {output}.")

    return df


def compare_mutation_lists(
    df1,
    df2,
    comp_cols=["sample", "Chr", "Start", "End"],
    keep_cols=[],
    keep_synonymous=False,
    exclude_samples=[],
):
    """
    takes 2 mutation dataframes and finds the overlap of mutations and marks these in the newly added column "comp" as exclusive
    args:
        comp_cols: columns to be used for defining equality
        keep_cols: columns from either df to be transfered to the other in case the rows match in comp_cols
        keep_synomymous: if False, rows with "ExonicFunc" == "synonymous"| "synonymous SNV" will be excluded; requires column "ExonicFunc"
        exclude_samples: list of samples to be exluded from comparison and output; requires column "sample"
    """

    # clean the files using args
    def filter_df(*dfs):
        """
        filter the list using the args
        """
        filter_dfs = []
        for df in dfs:
            if not keep_synonymous:
                df = df.loc[~df["ExonicFunc"].str.startswith("syn")]
            if exclude_samples:
                df = df.query("sample not in @exclude_samples")
            filter_dfs.append(df)
        return filter_dfs

    df1, df2 = filter_df(df1, df2)

    indicator_df = (
        df1.merge(df2, on=comp_cols, how="outer", indicator=True)
        .rename(dict(_merge="comp"), axis=1)
        .loc[:, comp_cols + ["comp"] + keep_cols]
    )

    df1 = df1.merge(indicator_df, how="left")
    df1["comp"] = df1["comp"].astype(str)
    df1.loc[df1["comp"] == "left_only", "comp"] = "exclusive"

    df2 = df2.merge(indicator_df, how="left")
    df2["comp"] = df2["comp"].astype(str)
    df2.loc[df2["comp"] == "right_only", "comp"] = "exclusive"

    return df1, df2


############# I/O ################################################


def load_file(file_name, chrom_sort=True):
    if file_name.endswith(".gz"):
        df = pd.read_csv(file_name, sep="\t", compression="gzip")
    else:
        df = pd.read_csv(file_name, sep="\t")

    if chrom_sort:
        df["Chr"] = pd.Categorical(
            df["Chr"], [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        )

    # sort by available columns
    df = df.sort_values(
        [col for col in ["sample", "Chr", "Start", "End"] if col in df.columns]
    )
    return df


def save_file(df, out_file, sort=False):

    if sort:
        df["Chr"] = pd.Categorical(
            df["Chr"], [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        )
        # sort by available columns
        df = df.sort_values(
            [col for col in ["sample", "Chr", "Start", "End"] if col in df.columns]
        )

        if out_file.endswith(".gz"):
            df.to_csv(out_file, sep="\t", index=False, compression="gzip")
        else:
            df.to_csv(out_file, sep="\t", index=False)

    return df


############## COMPUTATIONS #########################


def make_blocks(
    df,
    expand_to_previous=True,
    join_col="join",
    remove_join_col=True,
    block_col_name="block",
):
    """
    IMPORTANT!!!
    takes a boolean (or 1/0) join_col and returns blocks for contiguous rows with similar
    values for join_col
    arguments:
    expand_to_previous: if set to true, blocks of similar join_col will be expanded to the row below
    (I implement join criteria with shift, where only the second of two rows will be marked)
    """

    if expand_to_previous:
        # join to shift row
        block_start = df[join_col].astype(bool).shift(-1, fill_value=False) & ~(
            df[join_col].astype(bool)
        )
    else:
        # find first in block
        block_start = df[join_col].astype(bool) & ~(
            df[join_col].astype(bool).shift(1, fill_value=False)
        )
    # apply the blocks by summing the block_starts
    df[block_col_name] = block_start.astype(int).cumsum()
    # remove block numbers for non-block rows
    df[block_col_name] = df[block_col_name] * (df[join_col] | block_start)
    if remove_join_col:
        cols = [col for col in df.columns if not col == join_col]
        df = df.loc[:, cols]
    return df