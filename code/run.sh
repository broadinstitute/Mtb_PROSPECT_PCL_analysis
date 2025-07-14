#!/usr/bin/env bash
set -ex

# This is the master script for the capsule. When you click "Reproducible Run", the code in this file will execute.

matlab -nodisplay -r "ver; addpath(genpath('.'), '-begin'); subset_moas; moa_concordance_analysis; run_spectral_clustering; run_pcl_similarity_scoring; annotate_results; run_pcl_confidence_scoring; quit"

# Path to final RMarkdown file from original, reported cluster results found using Matlab 2020a
RMD_FILE="/code/combine_loocv_and_make_moa_predictions_using_original_pcl_clusters.Rmd"

# Final output file name and format
OUTPUT_FILE="/results/combine_loocv_and_make_moa_predictions_using_original_pcl_clusters_output.html"

# Run the final RMarkdown file
Rscript -e "rmarkdown::render('$RMD_FILE', output_file = '$OUTPUT_FILE', clean = TRUE)"

# Path to final RMarkdown file from demo cluster results found locally using Matlab 2020b
RMD_FILE="/code/combine_loocv_and_make_moa_predictions_using_demo_pcl_clusters.Rmd"

# Final output file name and format
OUTPUT_FILE="/results/combine_loocv_and_make_moa_predictions_using_demo_pcl_clusters_output.html"

# Run the final RMarkdown file
Rscript -e "rmarkdown::render('$RMD_FILE', output_file = '$OUTPUT_FILE', clean = TRUE)"