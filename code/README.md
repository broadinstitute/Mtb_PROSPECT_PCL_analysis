# MTb PROSPECT Perturbagen CLass (PCL) Analysis

MTb PROSPECT Perturbagen CLass (PCL) Analysis is a computational, reference-based method to predict MOA from Mycobacterium tuberculosis PROSPECT data by comparing chemical-genetic interaction profiles of novel compounds with those of a reference set of compounds. The software used when the code and data were originally built was Matlab 2020a with pre-installed toolboxes including Bioinformatics Toolbox, Parallel Computing Toolbox, and Statistics and Machine Learning Toolbox and R 4.1. The environment setup here to demo the analysis uses Matlab 2020b with pre-installed toolboxes and R 4.4. 

## Key notes: 

Matlab 2020b featured changes to the built-in sum, eigs function and others that slightly alter the cluster results from the Matlab spectralcluster function even when controlling for the random seed, number of workers for multi-threading, and order of treatments in the input matrices. The k-means++ clustering step of spectral clustering is highly sensitive to the randomized initialization of K clusters from the eigenvectors matrix especially for larger MOAs with greater number of K clusters. Controlling for the previously mentioned three factors leads to identical cluster (and downstream) results from our original findings when run using Matlab 2020a, but not Matlab 2020b where, due to changes in various functions, randomized initialization or the values of the eigenvector matrix marginally differ for some MOAs. Nonetheless, the cluster, PCL, and downstream results run here using Matlab 2020b differ only marginally from our findings when run originally using Matlab 2020a and repeat/reproduce with every run. Accordingly, cluster results from running with Matlab 2020a are included in the /data/ folder and used for all successive steps to show and recapitulate our reported findings. All steps are run and results are also provided for clusters resulting from processing here with Matlab 2020b. 

The leave-one-out cross-validation (LOOCV) over all 437 KABX compounds is memory and computationally intensive and was originally processed in parallel using multi-core processing and compute clusters which is unavailable on this platform. Here we demo the LOOCV results for one compound as to not exceed allotted compute time and storage and include in the /data/ folder the combined LOOCV results (as originally processed) for all KABX compounds for summary in the final Rmarkdown scripts.

## Abstract

PROSPECT (PRimary screening Of Strains to Prioritize Expanded Chemistry and Targets) is an antimicrobial discovery platform that measures chemical-genetic interactions between small molecules and a pool of bacterial mutants, each depleted of a different essential protein target, to identify active compounds while simultaneously providing insight into their mechanisms of action (MOA). We developed a computational, reference-based method termed Perturbagen CLass (PCL) analysis to predict MOA from Mycobacterium tuberculosis PROSPECT data by comparing chemical-genetic interaction profiles of novel compounds with those of a reference set of compounds. We assembled a diverse reference set of 437 compounds with known or suspected anti-tubercular activity and published MOA annotation. Using this reference set, PCL analysis had 74% sensitivity and 80% precision in leave-one-out cross-validation for active compounds whose MOA was well-represented in the reference set. When applied to a blinded set of 173 antitubercular leads from Glaxo Smith Kline, PCL analysis made high-confidence MOA predictions for 59 of the 71 compounds with previously assigned MOA based on published, experimental evidence, agreeing with these assignments 88% of the time. From the remaining 102 compounds without known MOAs, we made 61 novel MOA predictions including a range of cell wall targets and new, structurally diverse inhibitors of respiration. Finally, we identified a novel scaffold that lacked wild-type activity in the primary screen of an unbiased chemical library but, based on its PROSPECT chemical-genetic interaction profiles, PCL analysis predicted it to be a respiration, QcrB inhibitor. We confirmed the accuracy of this prediction while chemically optimizing this scaffold to achieve wild-type activity. Thus, PCL analysis of PROSPECT data can rapidly assign MOA to molecules previously identified with antimicrobial activity and to novel scaffolds from unbiased libraries that initially lack significant wild-type activity in a strategy to yield new, potent antitubercular compounds with annotated MOA.

## Running on Code Ocean

Please simply click the Reproducible Run button, then it will automatically run the PCL analysis and give the final predictions in results/fep/chembl.

## Key output/results for uploaded cluster results from Matlab 2020a run (original, reported findings):

/results/combine_loocv_and_make_moa_predictions_using_matlab_2020a_clusters_output.html
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/test_cmpd_pcl_based_moa_predictions_simplified_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/kabx_loocv_pcl_based_moa_predictions_simplified_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/kabx_demo_loocv_pcl_based_moa_predictions_simplified_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/pcls.gmt
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/pcl_cluster_members_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/pcl_cluster_members_summary_table.csv
spectral clustering summary .png plots for each MOA in /data/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/

## Key output/results for cluster results from local Matlab 2020b run:
/results/combine_loocv_and_make_moa_predictions_using_matlab_2020b_clusters_output.html
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/test_cmpd_pcl_based_moa_predictions_simplified_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/kabx_loocv_pcl_based_moa_predictions_simplified_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/kabx_demo_loocv_pcl_based_moa_predictions_simplified_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/pcls.gmt
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/pcl_cluster_members_table.csv
/results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/pcl_cluster_members_summary_table.csv
spectral clustering summary .png plots for each MOA in /results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/

## Code citations/dependencies
ConcensusGLM is available at http://github.com/eachanjohnson/concensusGLM
Eachan O. Johnson et. al., Large-scale chemical-genetics yields new M. tuberculosis inhibitor classes, Nature, July 2019, doi: 10.1038/s41586-019-1315-z
- not used in this capsule but was used in pre-processing of data to get L2FC for each strain x treatment which was later transformed and normalized to GR and sGR

CmapM is available on GitHub at https://github.com/cmap/cmapM
CmapR is available through Bioconductor and on GitHub at. https://github.com/cmap/cmapR
Subramanian, A. et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell (2017)
Enache, O. M. et al. The GCTx format and cmapPy, R, M, J packages: resources for optimized storage and integrated traversal of annotated dense matrices. Bioinformatics 35, 1427-1429 (2019). https://doi.org:10.1093/bioinformatics/bty784

## Usage (see run.sh)

set -ex

matlab -nodisplay -r "addpath(genpath('.'), '-begin'); subset_moas; moa_concordance_analysis; run_spectral_clustering; run_pcl_similarity_scoring; annotate_results; run_pcl_confidence_scoring; quit"

RMD_FILE="/code/combine_loocv_and_make_moa_predictions_using_matlab_2020a_clusters.Rmd"

OUTPUT_FILE="/results/combine_loocv_and_make_moa_predictions_using_matlab_2020a_clusters_output.html"

Rscript -e "rmarkdown::render('$RMD_FILE', output_file = '$OUTPUT_FILE', clean = TRUE)"

RMD_FILE="/code/combine_loocv_and_make_moa_predictions.Rmd"

OUTPUT_FILE="/results/combine_loocv_and_make_moa_predictions_using_matlab_2020b_clusters_output.html"

Rscript -e "rmarkdown::render('$RMD_FILE', output_file = '$OUTPUT_FILE', clean = TRUE)"

## Modifiable Inputs

Located as top sections of each of the Matlab scripts and RMarkdown files and should be carefully changed in each. 

### Defaults below:

#### General input and output filenames that should remain unchanged or if changed be reflected in each script accordingly

- not listed for length, see scripts

#### General LOOCV inputs (to be changed in each script accordingly)

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'number' % 'number' or 'list'

demo_loocv_number_cmpds = 1

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359','BRD-K87202646','BRD-K59853741', 'BRD-K27302037'} % Ciprofloxacin, Rifampin, Isoniazid, Q203, Thioacetazone

#### Spectral Clustering inputs (to be changed in each successive script accordingly)

rng_seed = 0; % specified seed for initializing the random number generator, Matlab factory default is the Mersenne Twister generator with seed 0 (see spectral_clustering_for_pcls.m for additional information)

num_threads_for_multithreading = 4 % number of computational threads to use for multi-threading enabled Matlab functions including eigs; data was originally processed with 4 CPU cores or threads, this is set for reproducibility across CPU hardware

thrsh_rank = 20 % threshold for average pairwise rank of correlation across KABX to connect treatments as mutual nearest-neighbors

thrsh_factor = 1; % factor to multiply by thrsh_rank if dynamic_thrsh_per_moa is false; default is 1

dynamic_thrsh_per_moa = false % if true then threshold is round(log(size of MOA) * thrsh_rank), otherwise identical threshold for every MOA

k_type = 'k_med_gap_den' % eigengap heuristic to take for estimating number of K clusters: k_num_zero, k_num_zero_plus_one, k_med_gap_den, k_gap_den (see create_laplacian_matrix.m for additional information)

show_hclust = true % if true then apply hierarchical clustering to correlation matrix prior to running spectral clustering, whether or not this is performed can minimally impact final cluster results especially for larger MOAs with greater number of K clusters due to stochasticity (randomized initialization) of k-means++ clustering; original results did apply hierarchical clustering but this step is not required for the method to be successful

#### Shared PCL similarity scoring, annotating results, and PCL confidence scoring inputs (to be changed in each accordingly)

% input whether to evaluate original cluster results produced using Matlab 2020a

use_matlab_2020a_clusters = true % Matlab 2020b featured changes to the built-in sum, eigs function and others that slightly alter the cluster results even when controlling for the random seed, number of workers for multi-threading, and order of treatments in the input matrices

matlab_2020a_clusters_outdir_name = 'clusters_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a'

matlab_2020a_clusters_pcls_results_outdir_name = 'pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a'

#### PCL similarity scoring inputs

min_clust_size = 2 % minimum cluster size to be considered for PCLs

#### PCL confidence scoring inputs

high_confidence_pcl_confidence_score_thres = 1; % defined as such for uniformity and simplicity across all PCLs

#### Combine LOOCV and Make MOA Predictions inputs (both Rmarkdown files for Matlab 2020a cluster results and local results)

any_dose_gr_filename = 'any_dose_min_gr.rds'

max_dose_gr_filename = 'max_dose_min_gr.rds'

- Filenames in /data/ folder with minimum curve-fit strain GR for each treatment (used to estimate compound activity)

## License

[MIT](https://choosealicense.com/licenses/mit/)