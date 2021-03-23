import pandas as pd
import os
from script_utils import show_output


def sort_filter2_df(
    df, cols={"Chr": True, "TVAF": False, "NVAF": False, "Start": True}
):
    """
    helper for sorting dfs for chromosomes using Chr, Start + cols in cols
    """
    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i}" for i in range(23)] + ["chrX", "chrY"]
    df["Chr"] = pd.Categorical(df["Chr"], chrom_list)
    return df.sort_values(list(cols.keys()), ascending=list(cols.values()))


def get_filter2(
    filter1_df,
    filter_settings,
    sample_df,
    keep_syn=False,
):

    # ########### LOADING FILTERS

    show_output(f"{len(filter1_df.index)} mutations found.", time=False)

    # remove syngeneic mutations if keep_syn is active (only valid for filter1)
    if keep_syn:
        is_exonic = ~filter1_df["ExonicFunc"].isin(["unknown", "synonymous SNV"])
        filter1_df = filter1_df[is_exonic]

    # use these population columns for checking PopFreq
    # could be refactored into params
    pop_cols = ["gnomAD_exome_ALL", "esp6500siv2_all", "dbSNP153_AltFreq"]

    output_base = filter2_output.replace(".loose.csv", "")

    # ######## FILTER 2 ######################
    def filter2(df, _filter="moderate"):

        # get thresholds
        show_output(f"filter: filter2-{_filter}")
        thresh = filter_settings.loc[f"filter2-{_filter}", :]
        # DEFINE CANDIDATE
        # used for rescue and special thresholds
        is_candidate = (
            (df["isCandidate"] == 1) | (df["isDriver"] == 1) | (df["ChipFreq"] > 0)
        )

        # #### SWITCH FOR AML7
        if "AML7" in filter_name:
            show_output("Including mutations on 7q to candidate list.", time=False)
            is7q = df["cytoBand"].str.contains("^7q")
            is_candidate = is_candidate | is7q

        # ##### TUMOR DEPTH ############
        tumor_depth = (df["TR2"] >= thresh["variantT"]) & (
            df["Tdepth"] >= thresh["Tdepth"]
        )

        # ##### VAF ##################
        # #### TVAF
        # either cand with lower TVAF or threshold TVAF
        TVAF = (is_candidate & (df["TVAF"] >= thresh["TVAF4Cand"])) | (
            df["TVAF"] >= thresh["TVAF"]
        )
        # ##### NVAF
        # NVAF is computed from upper threshold and a max proximity to TVAF (VAFSim)
        NVAF = (df["NVAF"] <= (thresh["VAFSim"] * df["TVAF"])) & (
            df["NVAF"] <= thresh["NVAF"]
        )

        # ##### EB/PoN-Filter ##########
        eb = (df["EBscore"] >= thresh["EBscore"]) if thresh["EBscore"] else True

        pon_eb = (eb & (df["PoN-Ratio"] < thresh["PoN-Ratio"])) | (
            df["PoN-Alt-NonZeros"] < thresh["PoN-Alt-NonZeros"]
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
            no_strand_polarity = no_strand_polarity = (
                df["TR2+"] / df["TR2"] <= pol
            ) & (df["TR2+"] / df["TR2"] >= (1 - pol))
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
        show_output(f"{stringency} {len(filter2_df.index)}", time=False)
        dropped_candidates_df = df[~filter_criteria & is_candidate]
        list_len = len(filter2_df.index)

        return (
            sort_filter2_df(filter2_df),
            sort_filter2_df(dropped_candidates_df),
            list_len,
        )

    # ################ OUTPUT #############################################################
    # #CSV ####

    filter2_dfs = {}
    dropped_dfs = {}
    df_lengths = {}
    for stringency in ["loose", "moderate", "strict"]:
        (
            filter2_dfs[stringency],
            dropped_dfs[stringency],
            df_lengths[stringency],
        ) = filter2(filter1_df, _filter=stringency)
        output_file = f"{output_base}.{stringency}.csv"
        show_output(
            f"Writing filter2.{stringency} ({df_lengths[stringency]}) to {output_file}"
        )
        filter2_dfs[stringency].to_csv(output_file, sep="\t", index=False)
    # write dropped files
    drop_file = f"{output_base}.dropped.csv"
    show_output(
        f"Writing {len(dropped_dfs['loose'].index)} muts to {drop_file}.", time=False
    )
    dropped_dfs["loose"].to_csv(drop_file, sep="\t", index=False)

    if excel_output:
        excel_file = f"{output_base}.xlsx"
        with pd.ExcelWriter(excel_file) as writer:
            # filter1
            filter1_df.to_excel(writer, sheet_name="filter1", index=False)
            for stringency in ["loose", "moderate", "strict"]:
                filter2_dfs[stringency].to_excel(
                    writer, sheet_name=stringency, index=False
                )

            show_output(f"Writing combined filters to excel file {excel_file}.")
            # write dropped files
            dropped_dfs["loose"].to_excel(writer, sheet_name="dropped", index=False)

    # create the filterbam_table for the selected stringency to be used by filterbam
    if filterbam_output:
        # if filterbam_stringency is one of loose, moderate, strict:

        # filterbam files should be sorted only for Chrom and position
        # they are input for filterbed
        sort_cols = {"Chr": True, "Start": True}
        if filterbam_stringency in filter2_dfs:
            sort_filter2_df(filter2_dfs[filterbam_stringency], cols=sort_cols).to_csv(
                filterbam_output, sep="\t", index=False
            )
        # if stringency is all or any other, use combination of loose and dropped for filterbam
        else:
            all_df = pd.concat(
                [filter2_dfs["loose"], dropped_dfs["loose"]]
            ).drop_duplicates()
            sort_filter2_df(all_df, cols=sort_cols).to_csv(
                filterbam_output, sep="\t", index=False
            )
