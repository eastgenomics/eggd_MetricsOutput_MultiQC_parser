#!/usr/bin/env python3.10

import argparse
import re
import pandas as pd
import numpy as np


def parse_metricsoutput_file(input_file):
    """
    Extract only relevant lines from MetricsOutput.tsv and store in dataframe.
    Transposes table so that metrics are columns and samples are rows.

    Parameters
    ----------
    input_file : str
        filepath to MetricsOutput.tsv output from eggd_tso500

    Returns
    ----------
    df : pd.DataFrame
        modified dataframe
    """
    try:
        df = pd.read_csv(input_file, sep='\t', header=None)
    except FileNotFoundError:
        print("Require input MetricsOutput.tsv file not found.")

    # get the indexes for key headings for MetricsOutput.tsv
    anaysis_status_index = df[df.iloc[:, 0].fillna('').str.contains("Analysis Status")].index[0]
    dna_lib_qc_index = df[df.iloc[:, 0].fillna('').str.contains("DNA Library QC Metrics")].index[0]
    dna_expanded_index = df[df.iloc[:, 0].fillna('').str.contains("DNA Expanded Metrics")].index[0]
    rna_lib_qc_index = df[df.iloc[:, 0].fillna('').str.contains("RNA Library QC Metrics")].index[0]
    rna_expanded_index = df[df.iloc[:, 0].fillna('').str.contains("RNA Expanded Metrics")].index[0]

    # First, get everything after "DNA Library QC Metrics" (to remove guideline columns adn clean)
    dna_rna_metrics_full = df.iloc[dna_lib_qc_index+1:]

    # drop the LSL and USL guidelines
    guidelines = ['LSL Guideline', 'USL Guideline']
    guideline_columns = dna_rna_metrics_full.applymap(lambda x: any(s in str(x) for s in guidelines)).any()
    dna_rna_metrics = dna_rna_metrics_full.loc[:, ~guideline_columns]

    # drop the DNA and RNA extended metrics
    dna_rna_metrics = dna_rna_metrics.drop(index=range(dna_expanded_index, rna_lib_qc_index))
    dna_rna_metrics = dna_rna_metrics.loc[:rna_expanded_index-1]

    # get indexes of lines that have matching regex (matches some string in square brackets)
    regex_index = dna_rna_metrics[dna_rna_metrics.iloc[:, 0].fillna('').str.match(r"^\[(.*)\]\s*$")].index.tolist()

    # drop those matching regex and empty lines
    dna_rna_metrics = dna_rna_metrics.drop(index=regex_index)
    dna_rna_metrics = dna_rna_metrics.dropna(how='all')

    # drop all duplicate Metric (UOM) rows with sample names
    dna_rna_metrics = dna_rna_metrics.drop_duplicates()

    # Now we have a clean df of relevant DNA and RNA metrics, get analysis status and concat
    # filter for just the analysis stats and drop empty lines
    analysis_stats = df.iloc[anaysis_status_index+1:dna_lib_qc_index-1]
    analysis_stats = analysis_stats.dropna(how= "all", axis=1)

    # set first cells as empty -- keeping it consistent to be able to concat later
    dna_rna_metrics.iloc[0, 0] = ""
    analysis_stats.iloc[0, 0] = ""

    # set first column as header for both dfs
    dna_rna_metrics.columns = dna_rna_metrics.iloc[0]
    dna_rna_metrics = dna_rna_metrics[1:]

    analysis_stats.columns = analysis_stats.iloc[0]
    analysis_stats = analysis_stats[1:]

    # combibe both dfs and set the index
    parsed_df = pd.concat([analysis_stats, dna_rna_metrics])
    parsed_df.set_index('', inplace=True)

    # transpose df
    transposed_df = parsed_df.transpose()

    return transposed_df


def edit_column_headers(edited_df):
    """
    Edits the column headers to be compatible with custom content feature of MultiQC.
    Rename MEDIAN_INSERT_SIZE (bp) for DNA to MEDIAN_INSERT_SIZE_DNA and 
    MEDIAN_INSERT_SIZE (Count) for RNA to MEDIAN_INSERT_SIZE_RNA.
    Remove the units of each metric.

    Parameters
    ----------
    edited_df : pd.DataFrame
        modified dataframe

    Returns
    ----------
    edited_df : pd.DataFrame
        modified dataframe
    """
    try:
        # differentiate MEDIAN_INSERT_SIZE between DNA and RNA
        # if a sample type is not present, there will be MEDIAN_INSERT_SIZE (NA) instead
        if "MEDIAN_INSERT_SIZE (bp)" in edited_df.columns:
            edited_df.columns = edited_df.columns.str.replace(
                "MEDIAN_INSERT_SIZE (bp)", "MEDIAN_INSERT_SIZE_DNA"
            )
        else:
            print("There are no DNA samples in this run.")

        if "MEDIAN_INSERT_SIZE (Count)" in edited_df.columns:
            edited_df.columns = edited_df.columns.str.replace(
                "MEDIAN_INSERT_SIZE (Count)", "MEDIAN_INSERT_SIZE_RNA"
            )
        else:
            print("There are no RNA samples in this run.")

        # removing the metric units
        edited_df.columns = edited_df.columns.str.split().str[0]

        return edited_df
    except:
        print("Error with editing column headers")


def add_contamination_bool(full_df):
    """
    Add CONTAMINATION_SUMMARY column, where False when 
    CONTAMINATION_SCORE > 3106 and CONTAMINATION_P_VALUE > 0.049, else True.

    Parameters
    ----------
    full_df : pd.DataFrame
        modified dataframe

    Returns
    ----------
    full_df : pd.DataFrame
        final modified dataframe
    """
    # replace NA string with null value
    full_df.replace("NA", pd.NA, inplace=True)

    # contamination metrics converted to float for numerical filtering
    try:
        full_df.loc[:, "CONTAMINATION_SCORE"] = full_df.loc[
            :, "CONTAMINATION_SCORE"
        ].astype(float)
        full_df.loc[:, "CONTAMINATION_P_VALUE"] = full_df.loc[
            :, "CONTAMINATION_P_VALUE"
        ].astype(float)
    except:
        print("Could not convert CONTAMINATION_SCORE and CONTAMINATION_P_VALUE to floats")
        raise ValueError

    # specify conditions
    nan_values = (
        full_df["CONTAMINATION_SCORE"].isna()
        | full_df["CONTAMINATION_P_VALUE"].isna()
    )
    filters = (full_df["CONTAMINATION_SCORE"] > 3106) & (
        full_df["CONTAMINATION_P_VALUE"] > 0.049
    )

    # initialise the new column
    full_df["CONTAMINATION_SUMMARY"] = None
    # for samples with contamination values, fill CONTAMINATION_SUMMARY with
    # the opposite of filters bool
    # if RNA sample, CONTAMINATION_SUMMARY = None
    # if DNA sample with contamination below threshold, CONTAMINATION_SUMMARY = True
    # if DNA sample with contamination above threshold, CONTAMINATION_SUMMARY = False
    try:
        full_df.loc[~nan_values, "CONTAMINATION_SUMMARY"] = ~filters[~nan_values]
    except:
        print("Could not specify CONTAMINATION_SUMMARY based on contamination thresholds.")

    # add header to index column
    full_df.index.name = "Sample"

    return full_df


def df_to_tsv(final_df, dna_output_filename, rna_output_filename):
    """
    Fills NaN values, checks for DNA and RNA samples and if present,
    saves the DNA and RNA samples in separate tsv files for input to MultiQC.

    Parameters
    ----------
    final_df : pd.DataFrame
        final modified dataframe

    Returns
    ----------
    MetricsOutput_MultiQC_DNA.tsv : file
        edited MetricsOutput.tsv file with only DNA samples for MultiQC input
    MetricsOutput_MultiQC_RNA.tsv : file
        edited MetricsOutput.tsv file with only RNA samples for MultiQC input
    """
    # fill pandas NaN values with NA
    final_df.fillna("NA", inplace=True)

    # split df into DNA and RNA samples
    for index in final_df.index:
        if "D" in index:
            dna_only_df = final_df.filter(like="D", axis=0)
            dna_only_df.to_csv(dna_output_filename, sep="\t")
        if "R" in index:
            rna_only_df = final_df.filter(like="R", axis=0)
            rna_only_df.to_csv(rna_output_filename, sep="\t")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="MetricsOutput.tsv editor for MultiQC input"
    )

    parser.add_argument(
        "tsv_input", type=str, help="filepath to MetricsOutput.tsv file"
    )
    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    # specify output filenames
    dna_output_filename = "MetricsOutput_MultiQC_DNA.tsv"
    rna_output_filename = "MetricsOutput_MultiQC_RNA.tsv"

    parsed_file = parse_metricsoutput_file(args.tsv_input)
    edited_df = edit_column_headers(parsed_file)
    final_df = add_contamination_bool(edited_df)
    df_to_tsv(final_df, dna_output_filename, rna_output_filename)


if __name__ == "__main__":
    main()
