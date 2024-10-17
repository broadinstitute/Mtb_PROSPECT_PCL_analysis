#!/bin/bash

#$ -l h_vmem=64G
#$ -l h_rt=3:00:00
#$ -j y
#$ -R y
#$ -cwd
#$ -m beas
#$ -M abond@broadinstitute.org

cd /idi/cgtb/code_prep_for_code_ocean/code

source $HOME/.my.bashrc
#matlab -nodisplay -r "addpath(genpath('.')); subset_moas; moa_concordance_analysis; run_spectral_clustering; quit"
#matlab -nodisplay -r "addpath(genpath('.')); run_pcl_similarity_scoring; annotate_results; quit"
matlab -nodisplay -r "addpath(genpath('.'), '-begin'); run_pcl_confidence_scoring; quit"

#matlab -nodisplay -r "addpath(genpath('.'), '-begin'); subset_moas; moa_concordance_analysis; run_spectral_clustering; run_pcl_similarity_scoring; annotate_results; run_pcl_confidence_scoring; quit"

# Path to final RMarkdown file
RMD_FILE="/code/combine_loocv_and_make_moa_predictions.Rmd"

# Final output file name and format
OUTPUT_FILE="/results/combine_loocv_and_make_moa_predictions_output.html"

# Run the final RMarkdown file
Rscript -e "rmarkdown::render('$RMD_FILE', output_file = '$OUTPUT_FILE', clean = TRUE)"

# Report resource consumption because it's not reported by default
echo "------------------------------"
qstat -j $JOB_ID | grep '^usage'
