#!/usr/bin/env python3.10

import argparse
import re
import pandas as pd
import numpy as np


# does this go here as global variables?
dna_output_filename = "MetricsOutput_MultiQC_DNA.tsv"
rna_output_filename = "MetricsOutput_MultiQC_RNA.tsv"


def parse_metricsoutput_file(input_file):
    """
    Extract only relevant lines from MetricsOutput.tsv and store in dataframe.
    This code partly taken from https://github.com/eastgenomics/multiqc_plugins/blob/develop/seglh_plugin/modules/tso500/tso500.py

    Parameters
    ----------
    input_file : file
        MetricsOutput.tsv output from eggd_tso500

    Returns
    ----------
    df : pd.DataFrame
        modified dataframe
    """
    with open(input_file, 'r', encoding='UTF-8') as file:
        group, sample_names, all_data = '', [], []
        for line in file:
            # match data block header
            m = re.match(r'^\[(.*)\]\s*$', line)
            if line.startswith('#'):
                # comment 
                continue
            elif len(line) == 0 or re.match(r'^\s+$',line):
                # empty line (reset)
                group, sample_names = '', []
                continue
            elif m:
                # is a group header name from matched Regexp (section name in square brackets => section in MultiQC report)
                group = m.group(1)
                #print(group)
            elif group:
                if group in ['Header']:
                    # global metrics/data
                    pass
                elif group in ['Analysis Status']:
                    contamination_value = line.rstrip().split('\t')[0:]
                    all_data.append(contamination_value)
                elif group.startswith('DNA Library') or group.startswith('RNA Library'):
                    # DNA data line (header or data)
                    if line.startswith("Metric "):
                        # is the header line, extract sample names from row
                        sample_names = line.rstrip().split('\t')[3:]
                        pass
                    else:
                        metric_values = line.rstrip().split('\t')[:1] + line.rstrip().split('\t')[3:]
                        all_data.append(metric_values)

    #print(all_data)

    df = pd.DataFrame(all_data)

    return df



# def transpose_table
def transpose_table(df):
    """
    Transposes table so that metrics are columns and samples are rows

    Parameters
    ----------
    df : pd.DataFrame
        modified dataframe

    Returns
    ----------
    transposed_df : pd.DataFrame
        transposed dataframe
    """
    transposed_df = df.transpose()

    # set sample names as index
    transposed_df.set_index(0, inplace=True)

    # set metrics as column headings
    transposed_df.columns = transposed_df.iloc[0]
    transposed_df = transposed_df[1:]

    return transposed_df



# def remove spaces from metrics, remove the parentheses, rename MEDIAN_INSERT_SIZE for DNA and RNA
def edit_column_headers(transposed_df):
    """
    Edits the column headers to be compatible with custom content feature of MultiQC.
    Rename MEDIAN_INSERT_SIZE (bp) for DNA to MEDIAN_INSERT_SIZE_DNA and MEDIAN_INSERT_SIZE (Count) 
    for RNA to MEDIAN_INSERT_SIZE_RNA.
    Remove the units of each metric.

    Parameters
    ----------
    transposed_df : pd.DataFrame
        modified dataframe
    
    Returns
    ----------
    transposed_df : pd.DataFrame
        modified dataframe
    """
    # differentiate MEDIAN_INSERT_SIZE between DNA and RNA
    transposed_df.columns = transposed_df.columns.str.replace('MEDIAN_INSERT_SIZE (bp)', 'MEDIAN_INSERT_SIZE_DNA')
    transposed_df.columns = transposed_df.columns.str.replace('MEDIAN_INSERT_SIZE (Count)', 'MEDIAN_INSERT_SIZE_RNA')
    
    # removing the metric units
    transposed_df.columns = transposed_df.columns.str.split().str[0]

    return transposed_df



# def add contamination true/false
def add_contamination_bool(transposed_df):
    """
    Add CONTAMINATION_SUMMARY column, where False when CONTAMINATION_SCORE > 3106 and 
    CONTAMINATION_P_VALUE > 0.049, else True.

    Parameters
    ----------
    transposed_df : pd.DataFrame
        modified dataframe
    
    Returns
    ----------
    result_df : pd.DataFrame
        final modified dataframe
    """
    # replace NA string with null value
    transposed_df.replace('NA', np.nan, inplace=True)

    # contamination metrics converted to float so filtering can be done
    transposed_df.loc[:, 'CONTAMINATION_SCORE'] = transposed_df.loc[:, 'CONTAMINATION_SCORE'].astype(float)
    transposed_df.loc[:, 'CONTAMINATION_P_VALUE'] = transposed_df.loc[:, 'CONTAMINATION_P_VALUE'].astype(float)

    # specify conditions
    nan_values = transposed_df['CONTAMINATION_SCORE'].isna() | transposed_df['CONTAMINATION_P_VALUE'].isna()
    filters = (transposed_df['CONTAMINATION_SCORE'] > 3106) & (transposed_df['CONTAMINATION_P_VALUE'] > 0.049)

    # initialise the new column
    transposed_df['CONTAMINATION_SUMMARY'] = None
    # for samples without NaN values, populate CONTAMINATION_SUMMARY with the opposite of filters bool
    transposed_df.loc[~nan_values, 'CONTAMINATION_SUMMARY'] = ~filters[~nan_values]

    # add header to index column
    transposed_df.index.name = "Sample"

    return transposed_df



# def save dataframe as tsv file
def df_to_tsv(result_df, dna_output_filename, rna_output_filename):
    """ # TODO: update this description
    Fills NaN values, separates DNA and RNA samples, and saves the modified dataframe 
    as a tsv file for input to MultiQC.

    Parameters
    ----------
    results_df : pd.DataFrame
        final modified dataframe

    Returns
    ----------
    MetricsOutput_MultiQC.tsv : file
        edited MetricsOutput.tsv file for MultiQC input
    """
    # fill pandas NaN values with NA
    result_df.fillna('NA', inplace=True)
    
    # split df into DNA and RNA samples
    for index in result_df.index:
        if "D" in index:
            dna_only_df = result_df.filter(like='D', axis=0)
            dna_only_df.to_csv(dna_output_filename, sep="\t")
        if "R" in index:
            rna_only_df = result_df.filter(like='R', axis=0)
            rna_only_df.to_csv(rna_output_filename, sep="\t")



def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(description='MetricsOutput.tsv \
                                     editor for MultiQC input')

    parser.add_argument('tsv_input', type=str,
                        help='filepath to MetricsOutput.tsv file')
    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    parsed_file = parse_metricsoutput_file(args.tsv_input)
    transposed_df = transpose_table(parsed_file)
    renamed_df = edit_column_headers(transposed_df)
    final_df = add_contamination_bool(renamed_df)
    df_to_tsv(final_df, dna_output_filename, rna_output_filename)


if __name__ == '__main__':
    main()
