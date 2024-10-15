% # Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, false);
%imatlab_export_fig('print-png')

% loocv inputs

prepare_loocv = false

demo_loocv = true

demo_loocv_number_or_list = 'number' % 'number' or 'list'

demo_loocv_number_cmpds = 2

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359','BRD-K87202646','BRD-K59853741', 'BRD-K27302037'} % Ciprofloxacin, Rifampin, Isoniazid, Q203, Thioacetazone

results_subdir_prefix = 'loocv_pcls/leave_out_cmpd_'

loocv_corr_filename = 'sGR_for_pcls_pearson_corr' % for faster runtime and lower storage, only will consider KABX dsCGI profiles in LOOCV

loocv_col_meta_all_filename = 'col_meta_kabx.txt' % for faster runtime and lower storage, only will consider KABX dsCGI profiles in LOOCV

unique_kabx_cmpds_tbl_path = '../results/kabx_pert_ids_tbl_for_loocv.txt'

% general inputs

clusters_gmt_filename = 'clusters_spectral_clust.gmt'

corr_filename = 'sGR_kabx_gsk_brd4310_pearson_corr'

col_meta_all_filename = 'col_meta.txt'

col_meta_kabx_filename = 'col_meta_kabx_for_pcls.txt'

save_out = false % save tabular file for each treatment's PCL similarity score

min_clust_size = 2

print_multi_target = false

stringify_cids = false

% previous Spectral Clustering inputs

thrsh_rank = 20 % threshold for average pairwise rank of correlation across KABX to connect treatments as mutual nearest-neighbors

dynamic_thrsh_per_moa = false % if true then threshold is round(log(size of MOA) * thrsh_rank), otherwise identical threshold for every MOA

k_type = 'k_med_gap_den' % eigengap heuristic to take for estimating number of K clusters: k_num_zero, k_num_zero_plus_one, k_med_gap_den, k_gap_den (see create_laplacian_matrix.m for additional information)

if dynamic_thrsh_per_moa
    prev_outdir_name = sprintf('clusters_spectral_clustering_thrsh_rank_le%dxlogsize_%s', thrsh_rank, k_type) % log(MOA size), i.e. the number of treatments/dsCGI profiles in the MOA
    outdir_name = sprintf('pcls_spectral_clustering_thrsh_rank_le%dxlogsize_%s', thrsh_rank, k_type) % log(MOA size), i.e. the number of treatments/dsCGI profiles in the MOA
else
    prev_outdir_name = sprintf('clusters_spectral_clustering_thrsh_rank_le%d_%s', thrsh_rank, k_type)
    outdir_name = sprintf('pcls_spectral_clustering_thrsh_rank_le%d_%s', thrsh_rank, k_type)
end

% # Run PCL similarity scoring 

clusters_path = fullfile(wkdir, prev_outdir_name, clusters_gmt_filename)
g = glob(fullfile(wkdir, [corr_filename,'_n*.gctx']));
c_path = g{1}
c_rank_path = []
col_meta_path = fullfile(wkdir, col_meta_all_filename)
col_meta_kabx_path = fullfile(wkdir, col_meta_kabx_filename)
outdir = fullfile(wkdir, outdir_name)

assert(exist(clusters_path) > 0)

assert(exist(c_path) > 0)

assert(exist(col_meta_path) > 0)

assert(exist(col_meta_kabx_path) > 0 )

pcl_similarity_scoring(clusters_path,c_path,c_rank_path,col_meta_path,col_meta_kabx_path,outdir,min_clust_size,print_multi_target,stringify_cids)

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
    
        clusters_path = fullfile(loo_wkdir, prev_outdir_name, clusters_gmt_filename)
        g = glob(fullfile(wkdir, [loocv_corr_filename,'_n*.gctx'])); % calculating PCL similarity scores for all KABX dsCGI profiles and excluding test compounds
        c_path = g{1}
        c_rank_path = []
        col_meta_path = fullfile(wkdir, loocv_col_meta_all_filename) % col_meta is for all of KABX only excluding test compounds
        col_meta_kabx_path = fullfile(loo_wkdir, col_meta_kabx_filename) % col_meta_kabx (in this iteration of LOOCV with a KABX compound treated as unknown) is the col_meta_kabx_for_pcls file in its subdirectory
        outdir = loo_outdir
        
        assert(exist(clusters_path) > 0)

        assert(exist(c_path) > 0)

        assert(exist(col_meta_path) > 0)

        assert(exist(col_meta_kabx_path) > 0 )
        
        pcl_similarity_scoring(clusters_path,c_path,c_rank_path,col_meta_path,col_meta_kabx_path,outdir,min_clust_size,print_multi_target,stringify_cids)
        
    end
    
end