# eggd_MetricsOutput_MultiQC_parser

## What does this app do?
This app produces tsv files that are formatted for input to MultiQC as tables. If both DNA and RNA samples are present, they are separated into 2 tsv files. This app uses MetricsOutput.tsv file as input. 

## What inputs are required for this app?

### Required
- `tsv_input` (`file`) - A single MetricsOutput.tsv report generated from the eggd_tso500 app

## How does this app work?
The app runs a python script to create two .tsv files from the provided MetricsOutput.tsv file. The output tsv files are modified for MultiQC to take as input in the form of tables. The general outline is as follows:

- Download the provided MetricsOutput.tsv file.
- Parses MetricsOutput.tsv to obtain only the samples and metrics information.
- Transpose the table so that metrics are columns and samples are rows.
- Edit column headers to remove spaces and units of each metric to be compatible with MultiQC custom content input. Also, differentiate between DNA and RNA MEDIAN_INSERT_SIZE.
- Add a contamination bool, where CONTAMINATION_SUMMARY = True when contamination metrics pass (below threshold), CONTAMINATION_SUMMARY = False when contamination metrics fail (above threshold), or CONTAMINATION_SUMMARY = None if RNA sample.
- Saves output tsv files, separated to DNA and RNA samples if both are present.


## What does this app output
Returns .tsv files with all the changes listed above. These .tsv files are then fed as input to MultiQC as tables, where it is called by the custom_content section in the TSO500 MultiQC config. 


## This app was created by East Genomics GLH
This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

