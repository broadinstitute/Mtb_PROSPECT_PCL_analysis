% # Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, false);
%imatlab_export_fig('print-png')

% loocv inputs

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'number' % 'number' or 'list'

demo_loocv_number_cmpds = 2

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359','BRD-K87202646','BRD-K59853741', 'BRD-K27302037'} % Ciprofloxacin, Rifampin, Isoniazid, Q203, Thioacetazone

results_subdir_prefix = 'loocv_pcls/leave_out_cmpd_'

loocv_save_fig = false % .png files of correlation and rank of correlation heatmaps for each MOA

loocv_save_gct = false % greater memory consumption than .gctx but can be opened in Morpheus, https://software.broadinstitute.org/morpheus

loocv_show_hclust = false % apply hierarchical clustering to MOA correlation and rank of correlation matrices before plotting

% inputs

unique_kabx_cmpds_tbl_path = '../results/kabx_pert_ids_tbl_for_loocv.txt'

corr_for_pcls_savepath = 'sGR_for_pcls_pearson_corr'

corr_rank_for_pcls_savepath = 'sGR_for_pcls_pearson_corr_rank'

moa_gmt_savepath = 'moas.gmt'

save_fig = true % .png files of correlation and rank of correlation heatmaps for each MOA

save_gct = false % greater memory consumption than .gctx but can be opened in Morpheus, https://software.broadinstitute.org/morpheus

show_hclust = false % apply hierarchical clustering to MOA correlation and rank of correlation matrices before plotting

% outputs

figdir = 'verify_moas/figures'

gctdir = 'verify_moas/gcts'

% # Run MOA concordance analysis - correlation and average rank of correlation across KABX for dsCGI profiles from each MOA

% # Subsets and saves individual gctx files for each MOA for quicker access during next step of spectral clustering and offers visualization of correlation and rank of correlation heatmaps

g = glob(fullfile(wkdir,strcat(corr_for_pcls_savepath,'_n*.gctx')));
grzs_corr = parse_gctx(g{1});

g = glob(fullfile(wkdir,strcat(corr_rank_for_pcls_savepath,'_n*.gctx')));
grzs_corr_rank = parse_gctx(g{1});

moas = parse_gmt(fullfile(wkdir, moa_gmt_savepath));

% function run_concordance_analysis(wkdir,c_corr,c_corr_rank,gmt,makefig,figdir,makegct,gctdir,show_hclust)
run_concordance_analysis(wkdir, grzs_corr, grzs_corr_rank, moas, save_fig, figdir, save_gct, gctdir, show_hclust)

disp('Done')

% # LOOCV section

if prepare_loocv

    unique_kabx_cmpds_tbl = rtable(unique_kabx_cmpds_tbl_path);

    size(unique_kabx_cmpds_tbl)
    headt(unique_kabx_cmpds_tbl)
    
    unique_kabx_cmpds_list = unique(unique_kabx_cmpds_tbl.kabx_cmpd);

    length(unique_kabx_cmpds_list)
    
    number_of_cmpds_loocv = length(unique_kabx_cmpds_list)
    
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

    for i = index_cmpds_loocv

        % If the current iteration number is a multiple of 50
        if mod(i, 50) == 0
            % Print a status message
            fprintf('Currently at iteration %d\n', i);
        end

        leave_out_cmpd = unique_kabx_cmpds_list(i);

        loo_wkdir = fullfile(wkdir, strcat(results_subdir_prefix, strjoin(unique_kabx_cmpds_list(i))));

        mk_cd_dir(loo_wkdir, false);
        
        g_loo = glob(fullfile(loo_wkdir,strcat(corr_for_pcls_savepath,'_n*.gctx')));
        grzs_corr_loo = parse_gctx(g_loo{1});

        g_loo = glob(fullfile(loo_wkdir,strcat(corr_rank_for_pcls_savepath,'_n*.gctx')));
        grzs_corr_rank_loo = parse_gctx(g_loo{1});
        
        moas_loo = parse_gmt(fullfile(loo_wkdir, moa_gmt_savepath));
        
        % function run_concordance_analysis(wkdir,c_corr,c_corr_rank,gmt,makefig,figdir,makegct,gctdir,show_hclust)
        run_concordance_analysis(loo_wkdir, grzs_corr_loo, grzs_corr_rank_loo, moas_loo, loocv_save_fig, figdir, loocv_save_gct, gctdir, loocv_show_hclust);
        
        disp('Done')
        
    end
end

