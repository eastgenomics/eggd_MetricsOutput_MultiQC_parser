#!/bin/bash
# eggd_MetricsOutput_MultiQC_parser
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.

# Exit at any point if there is any error and output each line as it is executed (for debugging)
set -e -x -o pipefail

main() {

    echo "Value of tsv_input: '$tsv_input'"

    # Download the provided MetricsOutput.tsv file
    dx download "$tsv_input" -o tsv_input

    # Install all python packages from this app
    python3 -m pip install --no-index --no-deps packages/*

    # Run python script with requires MetricsOutput.tsv file
    python3 edit_MetricsOutput.py tsv_input
    
    # Find output file from python script
    dna_tsv_file=$(find . -name "MetricsOutput_MultiQC_DNA.tsv")
    rna_tsv_file=$(find . -name "MetricsOutput_MultiQC_RNA.tsv")

    # Check and upload tsv files
    if [ -z "$dna_tsv_file" ]; then
        echo "There are no DNA samples in this run."
    else
        echo File "$dna_tsv_file" created.
        dna_output_file=$(dx upload "$dna_tsv_file" --brief)
        dx-jobutil-add-output dna_output_file "$dna_output_file" --class=file
    fi

    if [ -z "$rna_tsv_file" ]; then
        echo "There are no RNA samples in this run."
    else
        echo File "$rna_tsv_file" created.
        rna_output_file=$(dx upload "$rna_tsv_file" --brief)
        dx-jobutil-add-output rna_output_file "$rna_output_file" --class=file
    fi
}
