# Generate-CNV-report

## Introduction

Script to merge CNV and SVs from the Manta and ExomeDepth programs into a single file. Annotates with the gene the CNV/SV falls within.

## Requirements

python = 3.6
pandas = 0.23.4
scikit-allel = 1.2.0

## To Run

`python generateCNVReport.py --runid <run id> --output <output_file_name.csv>  --bed <custom_roi_bed> --exome_metrics <exome_depth_metrics_file> --manta_dir <manta_directory> --exome_dir <exome_depth_directory> --coverage_dir <coverage_dir>`

## Example

`python generateCNVReport.py --exome_metrics test/exome_depth/160722_M02641_0121_000000000-ARUB0_ExomeDepth_Metrics.txt --bed test/IlluminaTruSightCancer_CustomROI_b37.bed --manta_dir test/manta/ --coverage_dir test/depth/ --exome_dir test/exome_depth/  --output test.csv`

The expected output can be found within test/results.csv

