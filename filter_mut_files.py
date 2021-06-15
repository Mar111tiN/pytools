import os
import pandas as pd

# GLOBALS
chrom_list = [f"chr{i}" for i in range(23)] + ["chrX", "chrY"]


def sort_filter2_df(
    df, cols={"PatID": True, "Chr": True, "TVAF": False, "NVAF": False, "Start": True}
):
    """
    helper for sorting dfs for chromosomes using Chr, Start + cols in cols
    """
    # make Chr column categorical for sorting .. and sort
    df.loc[:, "Chr"] = pd.Categorical(df["Chr"], chrom_list)
    return df.sort_values(list(cols.keys()), ascending=list(cols.values()))


########### MUTATION FILTERS


def filter2(
    df,
    *,
    filter_settings,
    stringency="moderate",
    pop_cols=["gnomAD", "esp6500siv2_all", "dbSNP153_AltFreq"],
    sample_file="",
    NVAF_factor=0.25,
    isAML7=False,
):
    if sample_file:
        sample_df = pd.read_csv(sample_file, sep="\t")
        # apply normal_contamination to filter_df
        df = df.merge(
            sample_df.loc[:, ["PatID", "normal_contamination", "purity75"]]
        ).rename(dict(purity75="tumor_purity"), axis=1)
    else:
        df["normal_contamination"] = 0

    thresh = filter_settings.loc[f"filter2-{stringency}", :]
    # DEFINE CANDIDATE
    # used for rescue and special thresholds
    is_candidate = (
        (df["isCandidate"] == 1) | (df["AMLDriver"] == 1) | (df["ChipFreq"] > 0)
    )

    # #### SWITCH FOR AML7
    if isAML7:
        print("Including mutations on 7q to candidate list.")
        is7q = df["cytoBand"].str.contains("^7q")
        is_candidate = is_candidate | is7q

    # ##### TUMOR DEPTH ############
    tumor_depth = (df["TR2"] >= thresh["variantT"]) & (df["Tdepth"] >= thresh["Tdepth"])

    # ##### VAF ##################
    # #### TVAF
    # either cand with lower TVAF or threshold TVAF
    TVAF = (is_candidate & (df["TVAF"] >= thresh["TVAF4Cand"])) | (
        df["TVAF"] >= thresh["TVAF"]
    )
    # ##### NVAF
    # NVAF is computed from upper threshold and a max proximity to TVAF (VAFSim)
    NVAF = (df["NVAF"] <= (thresh["VAFSim"] * df["TVAF"])) & (
        df["NVAF"] <= thresh["NVAF"] + (df["normal_contamination"] * NVAF_factor)
    )

    # ##### EB/PoN-Filter ##########
    eb = (df["EBscore"] >= thresh["EBscore"]) if thresh["EBscore"] else True

    pon_eb = (eb & (df["PONRatio"] < thresh["PONRatio"])) | (
        df["PONAltNonZeros"] < thresh["PONAltNonZeros"]
    )

    # ############ HDR ####################
    HDR = (df["TumorHDRcount"] <= thresh["HDRcount"]) & (
        df["NormalHDRcount"] <= thresh["HDRcount"]
    )

    # ##### POPULATION #############
    if thresh["PopFreq"] == thresh["PopFreq"]:
        # init a noSNP series with True values for the looping
        noSNP = pd.Series(True, index=df.index)
        # go through pop_cols and filter..
        for col in pop_cols:
            # reformat population columns for filtering
            df.loc[df[col] == ".", col] = 0
            df[col] = df[col].fillna(0).astype(float)
            # combine the looped noSNP with the columns PopFreq checks
            noSNP = noSNP & (df[col] <= thresh["PopFreq"])
    else:
        noSNP = True

    # ####### STRANDBIAS / POLARITY ##########################
    # Strand Ratio (as FisherScore and simple)
    no_strand_bias = df["FisherScore"] <= thresh["FisherScore"]
    # Strand Polarity (filters out very uneven strand distribution of alt calls)
    if thresh.get("strandPolarity", None):
        pol = thresh["strandPolarity"]
        no_strand_polarity = no_strand_polarity = (df["TR2+"] / df["TR2"] <= pol) & (
            df["TR2+"] / df["TR2"] >= (1 - pol)
        )
    else:
        no_strand_polarity = True

    # let's check whether and is not too hard
    strandOK = no_strand_bias & no_strand_polarity

    # ########## RESCUE #####################
    # Clin_score is used for rescue of all mutations
    clin_score = df["ClinScore"] >= thresh["ClinScore"]
    rescue = clin_score

    rescue = False

    # ########### COMBINE CRITERIA ###############
    # rescue is only applied to disputable values within parens
    # criteria outside of parens are hard-filtered
    filter_criteria = (
        tumor_depth & pon_eb & NVAF & (noSNP & strandOK & TVAF & HDR | rescue)
    )

    filter2_df = df[filter_criteria]
    print(f"{stringency} {len(filter2_df.index)}")
    dropped_candidates_df = df[~filter_criteria & is_candidate]

    return (
        sort_filter2_df(
            filter2_df,
            cols={
                "PatID": True,
                "Chr": True,
                "TVAF": False,
                "NVAF": False,
                "Start": True,
            },
        ),
        sort_filter2_df(
            dropped_candidates_df,
            cols={
                "PatID": True,
                "Chr": True,
                "TVAF": False,
                "NVAF": False,
                "Start": True,
            },
        ),
    )


def filter2excel(filter2_df, filter2_dropped, out_file):
    """
    adds helpful additional sheets to the excel output of the filter file
    """

    # make additional agg files
    sample_df = (
        filter2_df.groupby(["PatID", "Gene"])
        .agg({"Chr": "count"})
        .reset_index("Gene")
        .groupby("PatID")
        .count()
        .loc[:, ["Gene"]]
        .reset_index()
    )

    gene_df = (
        filter2_df.groupby(["PatID", "Gene"])
        .first()
        .reset_index()
        .groupby("Gene")
        .agg({"PatID": "count"})
        .sort_values("PatID", ascending=False)
        .reset_index()
    )

    with pd.ExcelWriter(out_file) as writer:
        filter2_df.to_excel(writer, sheet_name="filter2", index=False)
        sample_df.to_excel(writer, sheet_name="samples", index=False)
        gene_df.to_excel(writer, sheet_name="genes", index=False)
        # print("dropped_df")
        # filter2_dropped.to_excel(writer, sheet_name="dropped", index=False)

    return filter2_df, filter2_dropped, sample_df


def write_filter(
    filter1_file,  # the filter1 file from the pipeline as csv
    filter_settings="",
    filter_name="AMLMono7",  # the sheet in the filter settings file
    stringencies=["loose", "moderate", "strict"],  # stringencies to be output
    sample_file="",  # file containing the sample information (see example_data/sample_list.txt)
    NVAF_factor=0.25,  # experimental factor for taking sample contamination into account
    csv_output=True,
    pop_cols=[],
    excel_output=True,
    skip_samples=[],
    isAML7=False,
):
    """
    takes a filter file (gz-compressed filter file)
    """
    if not (csv_output or excel_output):
        print("No written output")

    # load files
    filter_settings_df = pd.read_excel(
        filter_settings, sheet_name=filter_name, index_col=0, engine="openpyxl"
    )[:4]

    filter1_df = pd.read_csv(filter1_file, sep="\t", compression="gzip")
    if skip_samples:
        filter1_df = filter1_df.query("sample not in @skip_samples")

    filter_dfs = {"filter1": filter1_df}
    for stringency in stringencies:
        out_file_base = filter1_file.replace("filter1.csv.gz", "filter2")

        filter2_df, filter2_dropped = filter2(
            filter1_df,
            filter_settings=filter_settings_df,
            stringency=stringency,
            sample_file=sample_file,
            pop_cols=pop_cols,
            isAML7=isAML7,
        )
        if csv_output:
            filter2_df.to_csv(out_file_base + f".{stringency}.csv", sep="\t")
        if excel_output:
            _ = filter2excel(
                filter2_df, filter2_dropped, out_file_base + f".{stringency}.xlsx"
            )
        filter_dfs[stringency] = filter2_df
    return filter_dfs