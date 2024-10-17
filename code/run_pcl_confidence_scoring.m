% # Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, false);
%imatlab_export_fig('print-png')

% input whether to evaluate original cluster results produced using Matlab 2020a

use_matlab_2020a_clusters = false % Matlab 2020b featured changes to the built-in sum, eigs function and others that slightly alter the cluster results even when controlling for the random seed, number of workers for multi-threading, and order of treatments in the input matrices

matlab_2020a_clusters_outdir_name = 'clusters_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a'

matlab_2020a_clusters_pcls_results_outdir_name = 'pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a'

% loocv inputs

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'number' % 'number' or 'list'

demo_loocv_number_cmpds = 1

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359','BRD-K87202646','BRD-K59853741', 'BRD-K27302037'} % Ciprofloxacin, Rifampin, Isoniazid, Q203, Thioacetazone

results_subdir_prefix = 'loocv_pcls/leave_out_cmpd_'

unique_kabx_cmpds_tbl_path = '../results/kabx_pert_ids_tbl_for_loocv.txt'

% previous Spectral Clustering inputs

thrsh_rank = 20 % threshold for average pairwise rank of correlation across KABX to connect treatments as mutual nearest-neighbors

dynamic_thrsh_per_moa = false % if true then threshold is round(log(size of MOA) * thrsh_rank), otherwise identical threshold for every MOA

k_type = 'k_med_gap_den' % eigengap heuristic to take for estimating number of K clusters: k_num_zero, k_num_zero_plus_one, k_med_gap_den, k_gap_den (see create_laplacian_matrix.m for additional information)

if dynamic_thrsh_per_moa
    outdir_name = sprintf('pcls_spectral_clustering_thrsh_rank_le%dxlogsize_%s', thrsh_rank, k_type) % log(MOA size), i.e. the number of treatments/dsCGI profiles in the MOA
else
    outdir_name = sprintf('pcls_spectral_clustering_thrsh_rank_le%d_%s', thrsh_rank, k_type)
end

% step-specific inputs

pcls_filename = 'pcls.gmt'

col_meta_kabx_for_pcls_filename = 'col_meta_kabx_for_pcls.txt'

unknown_target_description_values = {'','NA','NAN','NaN','whole cell only','unknown'} % in case of different processing filetypes and NA values stored in import/export

make_fig = false

out_tbl_savename = 'by_pcl_similarity_to_confidence_score_thresholds_from_training_on_kabx.txt'

opt_tbl_savename = 'by_pcl_high_confidence_similarity_score_thresholds_from_training_on_kabx.txt'

out_tbl_test_cmpd_savename = 'by_pcl_similarity_to_confidence_score_test_cmpd_results.txt'

high_confidence_pcl_confidence_score_thres = 1; % defined as such for uniformity and simplicity across all PCLs

%pcl_stats_filename = 'figures_corr_summary.txt'

%pcl_stats_path = fullfile(wkdir, pcl_stats_filename)

% # PCL confidence scoring

col_meta_kabx_for_pcls_path = fullfile(wkdir, col_meta_kabx_for_pcls_filename)

exist(col_meta_kabx_for_pcls_path) > 0

assert(exist(col_meta_kabx_for_pcls_path) > 0, 'col_meta_kabx_for_pcls path does not exist')

% # Run PCL confidence scoring for uploaded Matlab 2020a clustering results

if use_matlab_2020a_clusters

    outdir = fullfile(wkdir, matlab_2020a_clusters_pcls_results_outdir_name)

    mk_cd_dir(outdir, false);

    pcls_path = fullfile(outdir, 'pcls.gmt')
    g = glob(fullfile(outdir,'pcl_similarity_score_n*.gctx'))
    similarity_score_path = g{1};
    g = glob(fullfile(outdir,'pcl_and_moa_agree_n*.gctx'))
    pcl_and_moa_agree_labels_path = g{1};
    outdir
    high_confidence_pcl_confidence_score_thres
    unknown_target_description_values
    make_fig
    out_tbl_savename
    opt_tbl_savename
    out_tbl_test_cmpd_savename

    exist(pcls_path) > 0

    exist(similarity_score_path) > 0

    exist(pcl_and_moa_agree_labels_path) > 0

    assert(exist(pcls_path) > 0, 'PCLs path does not exist')

    assert(exist(similarity_score_path) > 0, 'PCL similarity scores path does not exist')

    assert(exist(pcl_and_moa_agree_labels_path) > 0, 'pcl_and_moa_agree path does not exist')

    pcl_confidence_scoring(pcls_path,similarity_score_path,pcl_and_moa_agree_labels_path,col_meta_kabx_for_pcls_path,high_confidence_pcl_confidence_score_thres,outdir,unknown_target_description_values,make_fig,out_tbl_savename,opt_tbl_savename,out_tbl_test_cmpd_savename)

end

% # Run PCL confidence scoring for newly processed Matlab 2020b clustering results

outdir = fullfile(wkdir, outdir_name)

mk_cd_dir(outdir, false);

pcls_path = fullfile(outdir, 'pcls.gmt')
g = glob(fullfile(outdir,'pcl_similarity_score_n*.gctx'))
similarity_score_path = g{1};
g = glob(fullfile(outdir,'pcl_and_moa_agree_n*.gctx'))
pcl_and_moa_agree_labels_path = g{1};
outdir
high_confidence_pcl_confidence_score_thres
unknown_target_description_values
make_fig
out_tbl_savename
opt_tbl_savename
out_tbl_test_cmpd_savename

exist(pcls_path) > 0

exist(similarity_score_path) > 0

exist(pcl_and_moa_agree_labels_path) > 0

assert(exist(pcls_path) > 0, 'PCLs path does not exist')

assert(exist(similarity_score_path) > 0, 'PCL similarity scores path does not exist')

assert(exist(pcl_and_moa_agree_labels_path) > 0, 'pcl_and_moa_agree path does not exist')

pcl_confidence_scoring(pcls_path,similarity_score_path,pcl_and_moa_agree_labels_path,col_meta_kabx_for_pcls_path,high_confidence_pcl_confidence_score_thres,outdir,unknown_target_description_values,make_fig,out_tbl_savename,opt_tbl_savename,out_tbl_test_cmpd_savename)

% # LOOCV section

if prepare_loocv

    unique_kabx_cmpds_tbl = rtable(unique_kabx_cmpds_tbl_path);

    size(unique_kabx_cmpds_tbl)
    headt(unique_kabx_cmpds_tbl)
    
    unique_kabx_cmpds_list = unique(unique_kabx_cmpds_tbl.kabx_cmpd);

    length(unique_kabx_cmpds_list)
    
    number_of_cmpds_loocv = length(unique_kabx_cmpds_list)
    
    index_cmpds_loocv = 1:number_of_cmpds_loocv;
    
    if demo_loocv
       if strcmp(demo_loocv_number_or_list, 'number')
           number_of_cmpds_loocv = max(1, demo_loocv_number_cmpds)
           
           index_cmpds_loocv = 1:number_of_cmpds_loocv
           
       elseif strcmp(demo_loocv_number_or_list, 'list')
           number_of_cmpds_loocv = length(demo_loocv_list_cmpds)
           
           index_cmpds_loocv = find(ismember(unique_kabx_cmpds_list, demo_loocv_list_cmpds))'
       else
           error('Invalid input for demo_loocv_number_or_list: number or list')
       end
    end
    
    disp(sprintf('Number of KABX compounds to be processed in LOOCV: %d', length(index_cmpds_loocv)))
    
    for i = index_cmpds_loocv
    
        % If the current iteration number is a multiple of 50
        if mod(i, 50) == 0
            % Print a status message
            fprintf('Currently at iteration %d\n', i);
        end

        leave_out_cmpd = unique_kabx_cmpds_list(i);
        
        loo_wkdir = fullfile(wkdir, strcat(results_subdir_prefix, strjoin(unique_kabx_cmpds_list(i))));

        mk_cd_dir(loo_wkdir, false);
        
        loo_outdir = fullfile(loo_wkdir, outdir_name)
        
        mk_cd_dir(loo_outdir, false);
        
        % step specific commands
        
        pcls_path = fullfile(loo_outdir, 'pcls.gmt')
        g = glob(fullfile(loo_outdir,'pcl_similarity_score_n*.gctx'))
        similarity_score_path = g{1};
        g = glob(fullfile(loo_outdir,'pcl_and_moa_agree_n*.gctx'))
        pcl_and_moa_agree_labels_path = g{1};
        col_meta_kabx_for_pcls_path = fullfile(loo_wkdir, col_meta_kabx_for_pcls_filename)
        high_confidence_pcl_confidence_score_thres
        outdir = loo_outdir
        unknown_target_description_values
        make_fig
        out_tbl_savename
        opt_tbl_savename
        out_tbl_test_cmpd_savename
        
        exist(pcls_path) > 0

        exist(similarity_score_path) > 0

        exist(pcl_and_moa_agree_labels_path) > 0

        exist(col_meta_kabx_for_pcls_path) > 0

        assert(exist(pcls_path) > 0, 'PCLs path does not exist')

        assert(exist(similarity_score_path) > 0, 'PCL similarity scores path does not exist')

        assert(exist(pcl_and_moa_agree_labels_path) > 0, 'pcl_and_moa_agree path does not exist')

        assert(exist(col_meta_kabx_for_pcls_path) > 0, 'col_meta_kabx_for_pcls path does not exist')
        
        pcl_confidence_scoring(pcls_path,similarity_score_path,pcl_and_moa_agree_labels_path,col_meta_kabx_for_pcls_path,high_confidence_pcl_confidence_score_thres,outdir,unknown_target_description_values,make_fig,out_tbl_savename,opt_tbl_savename,out_tbl_test_cmpd_savename)
        
    end
    
end