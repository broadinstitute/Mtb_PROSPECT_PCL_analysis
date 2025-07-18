% ## Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, false);
%imatlab_export_fig('print-png')

% input whether to evaluate original cluster results produced using Matlab 2020a

use_matlab_2020a_clusters = true % if true process successive PCL analysis steps using original, spectral clustering cluster results from Matlab 2020a; Matlab 2020b introduced minor changes to built-in functions (e.g., sum, eigs) that can result in slight numerical differences in k-means++ initialization and spectral clustering outcomes. However, the overall impact on cluster assignments is minimal, and key downstream results remain highly consistent

matlab_2020a_clusters_outdir_name = 'original_clusters_spectral_clustering_thrsh_rank_le20_k_med_gap_den_matlab_2020a'

matlab_2020a_clusters_pcls_results_outdir_name = 'original_pcls_spectral_clustering_thrsh_rank_le20_k_med_gap_den'

% loocv inputs

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'list' % 'number' or 'list'

demo_loocv_number_cmpds = 1

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359'} % Ciprofloxacin, Rifampin

results_subdir_prefix = 'loocv_pcls/leave_out_cmpd_'

unique_reference_set_cmpds_tbl_path = '../results/reference_set_pert_ids_tbl_for_loocv.txt'

% general inputs

rr_path = fullfile(datadir, 'replicate_correlation.txt')

fgr_path = fullfile(datadir, 'fraction_gr.txt')

% previous Spectral Clustering inputs

thrsh_rank = 20 % threshold for average pairwise rank of correlation across reference set to connect conditions as mutual nearest-neighbors

dynamic_thrsh_per_moa = false % if true then threshold is round(log(size of MOA) * thrsh_rank), otherwise identical threshold for every MOA

k_type = 'k_med_gap_den' % eigengap heuristic to take for estimating number of K clusters: k_num_zero, k_num_zero_plus_one, k_med_gap_den, k_gap_den (see create_laplacian_matrix.m for additional information)

if dynamic_thrsh_per_moa
    outdir_name = sprintf('demo_pcls_spectral_clustering_thrsh_rank_le%dxlogsize_%s', thrsh_rank, k_type) % log(MOA size), i.e. the number of conditions/CGI profiles in the MOA
else
    outdir_name = sprintf('demo_pcls_spectral_clustering_thrsh_rank_le%d_%s', thrsh_rank, k_type)
end

exist(rr_path) > 0
exist(fgr_path) > 0

% ## Add replicate reproducibility 

rr = rtable(rr_path);

size(rr)
headt(rr)

% ## Add fraction of GR

fgr = rtable(fgr_path);

size(fgr)
headt(fgr)

rr.cid(~ismember(rr.cid, fgr.cid))

fgr.cid(~ismember(fgr.cid, rr.cid))

% # Process PCL cluster similarity scoring results from uploaded Matlab 2020a clusters

if use_matlab_2020a_clusters
    outdir = fullfile(wkdir, matlab_2020a_clusters_pcls_results_outdir_name)

    mk_cd_dir(outdir, false);

    % ## Full model section

    ls(outdir)

    % ### Add replicate reproducibility and fraction of GR to existing GCTs

    % Create gcts
    fields = {'cluster_median_corr','pcl_similarity_score','pcl_and_moa_agree'};

    for ii = 1:numel(fields)
        disp(fields{ii})
        g = glob(fullfile(outdir,strcat(fields{ii},'_n*.gctx')))
        s = parse_gctx(g{1});
        s.chd
        s.cid(~ismember(s.cid, rr.cid))
        s.cid(~ismember(s.cid, fgr.cid))
        s = annotate_ds(s,table2struct(rr(:,{'cid','rcorr','rcorr_rank'})),'dim','column','keyfield','cid');
        s = annotate_ds(s,table2struct(fgr(:,{'cid','x_median_total_count','frac_gr_le0_30'})),'dim','column','keyfield','cid');
        s.chd
        
        mkgctx(fullfile(outdir,strcat(fields{ii})), s)
    end

    s.chd
    s.rhd
end

% # Process PCL cluster similarity scoring results from newly processed Matlab 2020b clusters

outdir = fullfile(wkdir, outdir_name)

mk_cd_dir(outdir, false);

% ## Full model section

ls(outdir)

% ### Add replicate reproducibility and fraction of GR to existing GCTs

% Create gcts
fields = {'cluster_median_corr','pcl_similarity_score','pcl_and_moa_agree'};

for ii = 1:numel(fields)
    disp(fields{ii})
    g = glob(fullfile(outdir,strcat(fields{ii},'_n*.gctx')))
    s = parse_gctx(g{1});
    s.chd
    s.cid(~ismember(s.cid, rr.cid))
    s.cid(~ismember(s.cid, fgr.cid))
    s = annotate_ds(s,table2struct(rr(:,{'cid','rcorr','rcorr_rank'})),'dim','column','keyfield','cid');
    s = annotate_ds(s,table2struct(fgr(:,{'cid','x_median_total_count','frac_gr_le0_30'})),'dim','column','keyfield','cid');
    s.chd
    
    mkgctx(fullfile(outdir,strcat(fields{ii})), s)
end

s.chd
s.rhd

% ## LOOCV section

if prepare_loocv

    unique_reference_set_cmpds_tbl = rtable(unique_reference_set_cmpds_tbl_path);

    size(unique_reference_set_cmpds_tbl)
    headt(unique_reference_set_cmpds_tbl)
    
    unique_reference_set_cmpds_list = unique(unique_reference_set_cmpds_tbl.reference_set_cmpd);

    length(unique_reference_set_cmpds_list)
    
    number_of_cmpds_loocv = length(unique_reference_set_cmpds_list)
    
    index_cmpds_loocv = 1:number_of_cmpds_loocv;
    
    if demo_loocv
       if strcmp(demo_loocv_number_or_list, 'number')
           number_of_cmpds_loocv = max(1, demo_loocv_number_cmpds)
           
           index_cmpds_loocv = 1:number_of_cmpds_loocv
           
       elseif strcmp(demo_loocv_number_or_list, 'list')
           number_of_cmpds_loocv = length(demo_loocv_list_cmpds)
           
           index_cmpds_loocv = find(ismember(unique_reference_set_cmpds_list, demo_loocv_list_cmpds))'
       else
           error('Invalid input for demo_loocv_number_or_list: number or list')
       end
    end
    
    disp(sprintf('Number of reference set compounds to be processed in LOOCV: %d', length(index_cmpds_loocv)))
    
    for i = index_cmpds_loocv
    
        % If the current iteration number is a multiple of 50
        if mod(i, 50) == 0
            % Print a status message
            fprintf('Currently at iteration %d\n', i);
        end

        leave_out_cmpd = unique_reference_set_cmpds_list(i);
        
        loo_wkdir = fullfile(wkdir, strcat(results_subdir_prefix, strjoin(unique_reference_set_cmpds_list(i))));

        mk_cd_dir(loo_wkdir, false);
        
        loo_outdir = fullfile(loo_wkdir, outdir_name)
        
        mk_cd_dir(loo_outdir, false);
        
        % step specific commands
    
        ls(loo_outdir)
        
        % Add replicate reproducibility and fraction of GR to existing GCTs
        
        % Create gcts
        fields = {'cluster_median_corr','pcl_similarity_score','pcl_and_moa_agree'};

        for ii = 1:numel(fields)
            disp(fields{ii})
            g = glob(fullfile(loo_outdir,strcat(fields{ii},'_n*.gctx')))
            s = parse_gctx(g{1});
            s.chd
            s.cid(~ismember(s.cid, rr.cid))
            s.cid(~ismember(s.cid, fgr.cid))
            s = annotate_ds(s,table2struct(rr(:,{'cid','rcorr','rcorr_rank'})),'dim','column','keyfield','cid');
            s = annotate_ds(s,table2struct(fgr(:,{'cid','x_median_total_count','frac_gr_le0_30'})),'dim','column','keyfield','cid');
            s.chd

            mkgctx(fullfile(loo_outdir,strcat(fields{ii})), s)
        end
        
        s.chd
        s.rhd

    end
    
end