import os
from io import StringIO
from subprocess import PIPE, run
import pandas as pd
import numpy as np


def index_permutator(index_file):
    """
    takes an index.txt file containing all possible indices for a prep type and returns all possible combinations as an index_df
    """

    index_rename = {
        f"{i}{digit}{sep}{index}": f"i{digit}-Index"
        for i in "iI"
        for digit in [5, 7]
        for sep in "-_"
        for index in ["index", "Index"]
    }
    seq_rename = {
        f"{i}{digit}{sep}{seq}": f"i{digit}-Seq"
        for i in "iI"
        for digit in [5, 7]
        for sep in "-_"
        for seq in ["seq", "sequence", "Seq", "Sequence"]
    }
    df = (
        pd.read_csv(index_file, sep="\t")
        .rename(index_rename, axis=1)
        .rename(seq_rename, axis=1)
    )

    # auto-detect i5-columns
    hasI5 = False
    for col in df.columns:
        if "i5" in col or "I5" in col:
            hasI5 = True
            break
    index_df = df.iloc[:, [0, 1]].dropna()

    if hasI5:
        index_df["x"] = 0
        i5 = df.iloc[:, [2, 3]].dropna()
        i5["x"] = 0

        index_df = index_df.merge(i5).drop("x", axis=1)
        index_df["sample"] = index_df["i7-Index"] + "-" + index_df["i5-Index"]
        index_df = index_df.loc[
            :, ["sample", "i7-Index", "i7-Seq", "i5-Index", "i5-Seq"]
        ]
    else:
        index_df["sample"] = index_df["i7-Index"]
        index_df = index_df.loc[:, ["sample", "i7-Index", "i7-Seq"]]
    return index_df


def sheet2df(sample_sheet):
    """
    converts an Illumina sample sheet into an index_df
    """

    mawk_cmd = "$0~/\[Data\]/{data=1;next;}data"
    cmd = f"cat {sample_sheet} | mawk '{mawk_cmd}'"
    index_df = pd.read_csv(
        StringIO(run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode("utf-8")),
        sep=",",
    )
    index_df = index_df.rename(
        dict(
            Sample_ID="sample",
            I7_Index_ID="i7-Index",
            Index="i7-Seq",
            I5_Index_ID="i5-Index",
            Index2="i5-Seq",
            Sample_Project="Library",
        ),
        axis=1,
    )
    return index_df.loc[
        :, ["sample", "i7-Index", "i7-Seq", "i5-Index", "i5-Seq", "Library"]
    ]


def make_picard_sheets(
    index_df,
    library_name="NEBnext_IlluminaPrep",
    fastq_folder=".",
    output_file="",
    append_undetermined=True,
):
    """
    converts index_df into txt files for picard ExtractIlluminaBarcodes and IlluminaBarcodesToSam/Fastq
    <prefix>_barcode.txt is for ExtractIlluminaBarcodes
    <prefix>_multiplex.txt is for IlluminaBarcodesToSam/Fastq
    """

    # check if i5-columns are filled out
    hasI5 = "i5-Seq" in index_df.columns
    if hasI5:
        hasI5 = np.any((index_df["i5-Seq"] != "") & (~index_df["i5-Seq"].isna()))

    if hasI5:
        barcode_df = index_df.rename(
            {
                "i7-Seq": "barcode_sequence_1",
                "i5-Seq": "barcode_sequence_2",
                "sample": "barcode_name",
                "Library": "library_name",
            },
            axis=1,
        )

        # apply library name if not in list
        if not "library_name" in barcode_df.columns:
            barcode_df["library_name"] = library_name

        barcode_df = barcode_df.loc[
            :,
            [
                "barcode_sequence_1",
                "barcode_sequence_2",
                "barcode_name",
                "library_name",
            ],
        ]

        multiplex_df = barcode_df.rename(
            dict(
                barcode_name="OUTPUT_PREFIX",
                barcode_sequence_1="BARCODE_1",
                barcode_sequence_2="BARCODE_2",
            ),
            axis=1,
        ).loc[:, ["OUTPUT_PREFIX", "BARCODE_1", "BARCODE_2"]]

        # add the unplaces bam
        if append_undetermined:
            multiplex_df = multiplex_df.append(
                pd.Series(
                    {
                        "OUTPUT_PREFIX": "Undetermined",
                        "BARCODE_1": "N",
                        "BARCODE_2": "N",
                    }
                ),
                ignore_index=True,
            )

    else:  # no i5
        barcode_df = index_df.rename(
            {
                "i7-Seq": "barcode_sequence",
                "sample": "barcode_name",
                "Library": "library_name",
            },
            axis=1,
        )

        # apply library name if not in list
        if not "library_name" in barcode_df.columns:
            barcode_df["library_name"] = library_name

        barcode_df = barcode_df.loc[
            :, ["barcode_sequence", "barcode_name", "library_name"]
        ]

        multiplex_df = barcode_df.rename(
            dict(
                barcode_name="OUTPUT_PREFIX",
                barcode_sequence="BARCODE_1",
            ),
            axis=1,
        ).loc[:, ["OUTPUT_PREFIX", "BARCODE_1"]]
        # add the unplaced fastq
        if append_undetermined:
            multiplex_df = multiplex_df.append(
                pd.Series({"OUTPUT_PREFIX": "Undetermined", "BARCODE_1": "N"}),
                ignore_index=True,
            )

    multiplex_df["OUTPUT_PREFIX"] = (
        fastq_folder.rstrip("/") + "/" + multiplex_df["OUTPUT_PREFIX"]
    )
    # add the N column

    barcode_df.to_csv(f"{output_file}_barcode.txt", sep="\t", index=False)
    multiplex_df.to_csv(f"{output_file}_multiplex.txt", sep="\t", index=False)
    return barcode_df, multiplex_df
