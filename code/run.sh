#!/bin/bash

#$ -l h_vmem=64G
#$ -l h_rt=2:00:00
#$ -j y
#$ -R y
#$ -cwd
#$ -m beas
#$ -M abond@broadinstitute.org

cd /idi/cgtb/code_prep_for_code_ocean/code

source $HOME/.my.bashrc
#matlab -nodisplay -r "addpath(genpath('.')); subset_moas; moa_concordance_analysis; run_spectral_clustering; quit"
matlab -nodisplay -r "addpath(genpath('.')); run_pcl_similarity_scoring; annotate_results; quit"

# Report resource consumption because it's not reported by default
echo "------------------------------"
qstat -j $JOB_ID | grep '^usage'
