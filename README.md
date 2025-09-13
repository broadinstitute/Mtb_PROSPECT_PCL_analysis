# *Mycobacterium tuberculosis* PROSPECT Perturbagen CLass (PCL) Analysis

*Mycobacterium tuberculosis* (*Mtb*) PROSPECT **P**erturbagen **CL**ass (PCL) Analysis is a reference-based, computational method to predict mechanism-of-action (MOA) from *Mtb* PROSPECT data by comparing chemical-genetic interaction profiles of unknown compounds to those from a large reference set of compounds with known MOAs. The original analysis and results presented in our manuscript were generated using Matlab 2020a with the Bioinformatics Toolbox, Parallel Computing Toolbox, and Statistics and Machine Learning Toolbox, alongside R 4.1. The environment provided here (Code Ocean: https://doi.org/10.24433/CO.3013890.v1) for demonstrating the analysis runs Matlab 2020b-ubuntu20.04, the closest available version on Code Ocean, with the same pre-installed toolboxes and R 4.4.

## Abstract

We previously reported an antibiotic discovery screening platform that identifies whole-cell active compounds with high sensitivity while simultaneously providing mechanistic insight, necessary for hit prioritization. Named PROSPECT, (PRimary screening Of Strains to Prioritize Expanded Chemistry and Targets), this platform measures chemical-genetic interactions between small molecules and pooled Mycobacterium tuberculosis mutants, each depleted of a different essential protein. Here, we introduce Perturbagen CLass (PCL) analysis, a computational method that infers a compound's mechanism-of-action (MOA) by comparing its chemical-genetic interaction profile to those of a curated reference set of 437 known molecules. In leave-one-out cross-validation, we correctly predict MOA with 70% sensitivity and 75% precision, and achieve comparable results (69% sensitivity, 87% precision) with a test set of 75 antitubercular compounds with known MOA previously reported by GlaxoSmithKline (GSK). From 98 additional GSK antitubercular compounds with unknown MOA, we predict 60 to act via a reference MOA and functionally validate 29 compounds predicted to target respiration. Finally, from a set of ~5,000 compounds from larger unbiased libraries, we identify a novel QcrB-targeting scaffold that initially lacked wild-type activity, experimentally confirming this prediction while chemically optimizing this scaffold. PCL analysis of PROSPECT data enables rapid MOA assignment and hit prioritization, streamlining antimicrobial discovery.

## Acknowledgements

We thank Rob Bates and GlaxoSmithKline for kindly providing the TB set compounds for both primary screening and follow-up studies. RNA-Seq libraries were constructed and sequenced by the Infectious Disease and Microbiome Program’s Microbial Omics Core at the Broad Institute of MIT and Harvard. Funding for this work was provided by Bill and Melinda Gates Foundation (OPP1084233 D.T.H., INV-040933 D.T.H., INV-064678 D.T.H.), the Broad Institute Tuberculosis donor group, and the Pershing Square Foundation.

## System Requirements
- **Operating System:** Red Hat Enterprise Linux Server 7.9 (Maipo) (original analysis); Ubuntu 20.04 (local demo); other Linux-based systems may work
- **Matlab Versions:** Originally run on Matlab 2020a, demo on Matlab 2020b
- **Required Matlab Toolboxes:**  
    - Bioinformatics Toolbox  
    - Parallel Computing Toolbox  
    - Statistics and Machine Learning Toolbox  
- **R Versions:** Originally run on R 4.1, demo uses R 4.4
- **Hardware Requirements:**  
    - **Memory:** At least **20GB RAM** required per job  
    - **CPU:** Multi-core processor recommended
    - **Cluster Usage:** The scripts are optimized for running on a **high-performance computing (HPC) cluster** using **UGE/SGE**  
    - **For full LOOCV runs:** Parellelizing LOOCV iterations for each compound across multiple HPC cluster jobs is required for efficient processing
- **Expected Runtime:** Minimally **3 hours per job** if no LOOCV iterations are processed

## Installation Guide
- **Code Ocean users:** Simply click "Reproducible Run" to execute the pipeline.
- **Local installation (if applicable):**  
    - Install Matlab (2020a preferred for reproducing results, 2020b for demo).
    - Ensure required toolboxes are installed.
    - Install R (4.1 or 4.4) and required R dependencies (see packages in Environment).
    - Run postInstall script to download and initalize CmapM Matlab library
- **Typical install time:** ~1 hr on a standard desktop.

## Key notes: 

- The core computational pipeline remains robust and reproducible across Matlab versions.
- The clustering step in PCL analysis employs spectral clustering with k-means++ initialization, which is inherently sensitive to randomized cluster initialization, particularly for larger MOAs with a greater number of clusters.
- Matlab 2020b introduced minor changes to built-in functions such as sum and eigs, which can lead to small numerical differences in the eigenvector matrix. While this does not affect the overall performance of PCL analysis, it can result in minor variations in specific cluster assignments depending on the Matlab software version due to the sensitivity of k-means++ initialization.
- To ensure consistent clustering results within a given Matlab version, we controlled for key factors including random seed initialization, multi-threading parameters, and input data ordering. As a result, spectral clustering yields identical clusters when rerun within the same Matlab version.
- Critically, all downstream analyses, including MOA predictions and cross-validation results, remain stable and consistent across versions, demonstrating the robustness of the approach.
- To ensure direct comparability with our originally reported results, we provide the originally computed cluster assignments (from Matlab 2020a) in the /data/ folder, which are used for all downstream steps to recapitulate reported findings.
- Additionally, the pipeline runs using Matlab 2020b-generated clusters, and results remain highly consistent with the original analysis, further reinforcing the robustness of the methodology.

## Demonstration of Leave-One-Out Cross-Validation (LOOCV)

- Full LOOCV analysis across all 437 reference set compounds is computationally intensive and was originally executed in parallel using multi-core processing on high-performance compute clusters.
- Due to the resource constraints of this demo environment, we provide example LOOCV for two compounds, ciprofloxacin (BRD-K04804440) and rifampin (BRD-K01507359) and include precomputed LOOCV results (originally generated in Matlab 2020a) for all reference compounds in the /data/ folder.
- These precomputed results are summarized in the final Rmarkdown scripts, ensuring full transparency and reproducibility of the methodology.

## Running on Code Ocean

Simply click the Reproducible Run button to automatically execute the PCL analysis pipeline. Below are the key output files and a comparison of results generated using Matlab 2020a (original, reported findings) and those computed here as a demo using Matlab 2020b.

## Full Reproducibility
- **Original reported findings (Matlab 2020a)** are provided in `/results/original_*` for strict comparability.
- **Demo results (Matlab 2020b)** confirm the robustness of the method with **minimal differences**.
- Full LOOCV results **from the original analysis** are included in `/data/` for reference.

##  Key Output Directories

We use the following naming conventions to clearly differentiate between the original reported findings and the demo run for reproducibility:

- /results/original_* → Original results (using Matlab 2020a, as reported in our manuscript).
- /results/demo_* → Demo results (using Matlab 2020b, run here for reproducibility).

## Key Output Files

### Original Reported Findings (spectral clustering results using Matlab 2020a)

/results/**combine_loocv_and_make_moa_predictions_using_original_pcl_clusters_output.html**
- Includes accuracy statistics and plots on reference set in LOOCV and in full model (fit to data/training) at various data levels 
    - pert_id, unique identifier for compound based on structure
    - broad_id, unique identifier for commercial or synthetic lot of a compound
    - proj_broad_id, unique identifier for an instance of a specific lot of compound screened in a particular screening wave
    - cid/condition, unique screening wave instance of a specific lot of compound screened at a specific dose - CGI, chemical-genetic interaction

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**test_cmpd_pcl_based_moa_predictions_simplified_table.csv**
- GSK TB set and BRD4310 PCL-based MOA predictions from original cluster results processed using Matlab 2020a

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**reference_set_loocv_pcl_based_moa_predictions_simplified_table.csv**
- reference set PCL-based MOA predictions in LOOCV as processed originally using Matlab 2020a
- identical to file in /results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**reference_set_demo_loocv_pcl_based_moa_predictions_simplified_table.csv**
- PCL-based MOA predictions in LOOCV for subset of compounds in demo and processed locally using Matlab 2020b
- identical to file in /results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**pcls.gmt**
- GMT file with PCL cluster membership (CGI profiles in each PCL cluster from original spectral clustering using Matlab 2020a)

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**pcl_cluster_members_table.csv**
- PCL cluster membership table (from original spectral clustering using Matlab 2020a)

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**pcl_cluster_members_summary_table.csv**
- PCL cluster membership summary table (from original spectral clustering using Matlab 2020a)

spectral clustering summary .png plots for each MOA in /**data**/**original_clusters_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a**/
- for each MOA visualization of each of the steps of spectral clustering (as processed originally using Matlab 2020a)

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**by_pcl_similarity_to_confidence_score_thresholds_from_training_on_reference_set.txt.gz**
- compressed tab-separated file with full model PCL similarity score to confidence score mapping for each PCL cluster learned from reference set (from original spectral clustering using Matlab 2020a)

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**by_pcl_high_confidence_similarity_score_thresholds_from_training_on_reference_set.txt.gz**
- compressed tab-separated file with the high-confidence PCL similarity score thresholds for each PCL cluster learned from reference set (from original spectral clustering using Matlab 2020a)

/results/original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**by_pcl_similarity_to_confidence_score_test_cmpd_results.txt.gz**
- compressed tab-separated file with the PCL similarity score and estimated PCL confidence scores for every test compound CGI profile (condition) to each PCL cluster (from original spectral clustering using Matlab 2020a)
- these results are used in combine_loocv_and_make_moa_predictions_using_matlab_2020a.Rmd to make one MOA prediction for every test compound via a majority rules tiebreaker in which the PCL cluster that had the highest confidence score over the greatest number of doses (CGI profiles) and highest similarity score overall, if necessary, was selected. The predicted MOA was considered a high-confidence MOA assignment if the PCL confidence score was equal to 1. A compound’s MOA was considered “uncertain” if its highest PCL confidence score was below 1

### Demo Run (spectral clustering results using Matlab 2020b)

/results/**combine_loocv_and_make_moa_predictions_using_demo_pcl_clusters_output.html**
- Includes accuracy statistics and plots on reference set in LOOCV and in full model (fit to data/training) at various data levels 
    - pert_id, unique identifier for compound based on structure
    - broad_id, unique identifier for commercial or synthetic lot of a compound
    - proj_broad_id, unique identifier for an instance of a specific lot of compound screened in a particular screening wave
    - cid/condition, unique screening wave instance of a specific lot of compound screened at a specific dose - CGI, chemical-genetic interaction

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**test_cmpd_pcl_based_moa_predictions_simplified_table.csv**
- GSK TB set and BRD4310 PCL-based MOA predictions from cluster results processed locally using Matlab 2020b

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**reference_set_loocv_pcl_based_moa_predictions_simplified_table.csv**
- reference set PCL-based MOA predictions in LOOCV as processed originally using Matlab 2020a
- identical to file in /results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**reference_set_demo_loocv_pcl_based_moa_predictions_simplified_table.csv**
- PCL-based MOA predictions in LOOCV for subset of compounds in demo and processed locally using Matlab 2020b
- identical to file in /results/pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a/

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**pcls.gmt**
- GMT file with PCL cluster membership (CGI profiles in each PCL cluster from spectral clustering locally using Matlab 2020b)

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**pcl_cluster_members_table.csv**
- PCL cluster membership table (from spectral clustering locally using Matlab 2020b)

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**pcl_cluster_members_summary_table.csv**
- PCL cluster membership summary table (from spectral clustering locally using Matlab 2020b)

spectral clustering summary .png plots for each MOA in /results/**demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den**/
- for each MOA visualization of each of the steps of spectral clustering (as locally using Matlab 2020b)

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**by_pcl_similarity_to_confidence_score_thresholds_from_training_on_reference_set.txt.gz**
- compressed tab-separated file with full model PCL similarity score to confidence score mapping for each PCL cluster learned from reference set (from spectral clustering locally using Matlab 2020b)

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**by_pcl_high_confidence_similarity_score_thresholds_from_training_on_reference_set.txt.gz**
- compressed tab-separated file with the high-confidence PCL similarity score thresholds for each PCL cluster learned from reference set (from spectral clustering locally using Matlab 2020b)

/results/demo_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den/**by_pcl_similarity_to_confidence_score_test_cmpd_results.txt.gz**
- compressed tab-separated file with the PCL similarity score and estimated PCL confidence scores for every test compound CGI profile (condition) to each PCL cluster (from spectral clustering locally using Matlab 2020b)
- these results are used in combine_loocv_and_make_moa_predictions_using_matlab_2020a.Rmd to make one MOA prediction for every test compound via a majority rules tiebreaker in which the PCL cluster that had the highest confidence score over the greatest number of doses (CGI profiles) and highest similarity score overall, if necessary, was selected. The predicted MOA was considered a high-confidence MOA assignment if the PCL confidence score was equal to 1. A compound’s MOA was considered “uncertain” if its highest PCL confidence score was below 1

## Comparison of Results:

- Minimal differences between Matlab 2020a and Matlab 2020b spectral clustering
    - PCL-based MOA predictions remain highly consistent, demonstrating the robustness of the pipeline.
    - Matlab 2020b spectral clustering identified 1,139 PCL clusters vs. 1,140 from Matlab 2020a.
    - Accuracy statistics for reference set LOOCV and annotated GSK test set predictions differ only marginally.

- Why do slight differences occur?
    - The k-means++ initialization step in spectral clustering is inherently sensitive to variations in the input eigenvector matrix, which can arise from numerical precision differences or small perturbations in data ordering, particularly for larger MOAs with more k clusters. Minor numerical differences introduced in Matlab 2020b functions slightly alter a small number of spectral clustering outcomes across Matlab versions but do not impact overall MOA predictions or conclusions.
    - Within a given Matlab version, spectral clustering results are fully reproducible — we control for random seed initialization, multi-threading settings, and input data ordering, ensuring identical cluster assignments when rerun.

- Why does this not impact conclusions?
    - PCL-based MOA predictions remain stable, with over 97% agreement in top predicted PCLs across Matlab versions.
    - Key findings, including the identification of QcrB inhibitors, annotated GSK test set MOA prediction accuracy, and LOOCV performance, remain unchanged.
    - Precomputed Matlab 2020a clustering results are provided for strict reproducibility with the original reported findings, while Matlab 2020b results confirm consistency.

BRD4310 PCL-based MOA prediction
- no difference; high-confidence QcrB

Of 75 GSK with known, annotated MOA:
- 60 with high-confidence PCL-based MOA prediction (52 of 60 assigned MOAs concordant with annotated MOA = 87% precision; 52 / 75 = 69% sensitivity), 15 uncertain compared to 62 with high-confidence PCL-based MOA prediction (54 of 62 assigned MOAs concordant with annotated MOA = 87% precision; 54 / 75 = 72% sensitivity), 13 uncertain 

Of 98 GSK with unknown, unannotated MOA:
- 60 with high-confidence PCL-based MOA prediction, 38 uncertain compared to 59 with high-confidence PCL-based MOA prediction, 39 uncertain

## Code citations/dependencies
ConcensusGLM is available at http://github.com/eachanjohnson/concensusGLM
- not used in this capsule but was used in pre-processing of data to get L2FC for each strain x condition which was later transformed and normalized to GR and sGR

Eachan O. Johnson et. al., Large-scale chemical-genetics yields new M. tuberculosis inhibitor classes, Nature, July 2019, doi: 10.1038/s41586-019-1315-z

CmapM is available on GitHub at https://github.com/cmap/cmapM

CmapR is available through Bioconductor and on GitHub at https://github.com/cmap/cmapR

Subramanian, A. et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell (2017)

Enache, O. M. et al. The GCTx format and cmapPy, R, M, J packages: resources for optimized storage and integrated traversal of annotated dense matrices. Bioinformatics 35, 1427-1429 (2019). https://doi.org:10.1093/bioinformatics/bty784

## Usage (see run.sh)

```bash
#!/usr/bin/env bash
set -ex

matlab -nodisplay -r "ver; addpath(genpath('.'), '-begin'); subset_moas; moa_concordance_analysis; run_spectral_clustering; run_pcl_similarity_scoring; annotate_results; run_pcl_confidence_scoring; quit"

RMD_FILE="/code/combine_loocv_and_make_moa_predictions_using_original_pcl_clusters.Rmd"

OUTPUT_FILE="/results/combine_loocv_and_make_moa_predictions_using_original_pcl_clusters_output.html"

Rscript -e "rmarkdown::render('$RMD_FILE', output_file = '$OUTPUT_FILE', clean = TRUE)"

RMD_FILE="/code/combine_loocv_and_make_moa_predictions_using_demo_pcl_clusters.Rmd"

OUTPUT_FILE="/results/combine_loocv_and_make_moa_predictions_using_demo_pcl_clusters_output.html"

Rscript -e "rmarkdown::render('$RMD_FILE', output_file = '$OUTPUT_FILE', clean = TRUE)"
```

## Modifiable Inputs

Located as top sections of each of the Matlab scripts and RMarkdown files and should be carefully changed in each. 

### Defaults below:

#### General input and output filenames that should remain unchanged or if changed be reflected in each script accordingly

- not listed for length, see scripts

#### General LOOCV inputs (to be changed in each script accordingly)

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'list' % 'number' or 'list'

demo_loocv_number_cmpds = 1

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359'} % Ciprofloxacin, Rifampin

#### Spectral Clustering inputs (to be changed in each successive script accordingly)

rng_seed = 0; % specified seed for initializing the random number generator, Matlab factory default is the Mersenne Twister generator with seed 0 (see spectral_clustering_for_pcls.m for additional information)

num_threads_for_multithreading = 4 % number of computational threads to use for multi-threading enabled Matlab functions including eigs; data was originally processed with 4 CPU cores or threads, this is set for reproducibility across CPU hardware

thrsh_rank = 20 % threshold for average pairwise rank of correlation across reference set to connect conditions as mutual nearest-neighbors; default is 20

thrsh_factor = 1; % factor to multiply by thrsh_rank if dynamic_thrsh_per_moa is false; default is 1

dynamic_thrsh_per_moa = false % if true then threshold is round(log(size of MOA) * thrsh_rank), otherwise identical threshold for every MOA; was tested limitedly and found to achieve slightly worse LOOCV MOA prediction performance, default is false

k_type = 'k_med_gap_den' % eigengap heuristic to take for estimating number of K clusters: k_num_zero, k_num_zero_plus_one, k_med_gap_den, k_gap_den (see create_laplacian_matrix.m for additional information); default is k_med_gap_den

show_hclust = true % if true hierarchical clustering is applied to the correlation matrix before spectral clustering. While this step was used in the original analysis, it is not required for the method to be successful. For larger MOAs with more K clusters, hierarchical clustering may slightly influence final cluster assignments due to the inherent stochasticity of k-means++ initialization, but overall, its impact on results is minimal

#### Shared PCL similarity scoring, annotating results, and PCL confidence scoring inputs (to be changed in each accordingly)

% input whether to evaluate original cluster results produced using Matlab 2020a

use_matlab_2020a_clusters = true % if true process successive PCL analysis steps using original, spectral clustering cluster results from Matlab 2020a; Matlab 2020b introduced minor changes to built-in functions (e.g., sum, eigs) that can result in slight numerical differences in k-means++ initialization and spectral clustering outcomes. However, the overall impact on cluster assignments is minimal, and key downstream results remain highly consistent

matlab_2020a_clusters_outdir_name = 'original_clusters_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a'

matlab_2020a_clusters_pcls_results_outdir_name = 'original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den'

#### PCL similarity scoring inputs

min_clust_size = 2 % minimum cluster size to be considered for PCL clusters

#### PCL confidence scoring inputs

high_confidence_pcl_confidence_score_thres = 1; % defined as such for uniformity and simplicity across all PCL clusters

#### Combine LOOCV and Make MOA Predictions inputs (both Rmarkdown files for Matlab 2020a cluster results and local results)

loocv_opt_tbl_combined_full_shared_filename = "loocv_combined_by_pcl_high_confidence_similarity_score_thresholds_from_training_on_reference_set_tbls_full_shared.txt.gz"

loocv_test_cmpd_results_combined_full_shared_filename = "loocv_combined_test_cmpd_results_full_shared.txt.gz"

- Filenames in /data/ folder for precomputed LOOCV results over all reference set compounds

any_dose_gr_filename = 'any_dose_min_gr.txt'

max_dose_gr_filename = 'max_dose_min_gr.txt'

- Filenames in /data/ folder with minimum curve-fit strain GR for each condition (used to estimate compound activity)

## License

[MIT](https://choosealicense.com/licenses/mit/)
