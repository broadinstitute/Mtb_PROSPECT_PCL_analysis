% # Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, true);
%imatlab_export_fig('print-png')

% loocv inputs

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'list' % 'number' or 'list'

demo_loocv_number_cmpds = 1

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359'} % Ciprofloxacin, Rifampin

results_subdir_prefix = 'loocv_pcls/leave_out_cmpd_'

% inputs

reference_set_annot_path = '../data/reference_set_annotations.xlsx'

GR_filename = 'GR_reference_set_gsk_brd4310_n10819x340.gctx'

sGR_filename = 'sGR_reference_set_gsk_brd4310_n10819x340.gctx'

filter_moa_annotations = {'NaN','unknown','whole cell only'}

load_and_save_GR = false

save_gct = false % greater memory consumption than .gctx but can be opened in Morpheus, https://software.broadinstitute.org/morpheus

% outputs

col_meta_savepath = 'col_meta.txt'

col_meta_reference_set_savepath = 'col_meta_reference_set.txt'

col_meta_reference_set_for_pcls_savepath = 'col_meta_reference_set_for_pcls.txt'

corr_savepath = 'sGR_reference_set_gsk_brd4310_pearson_corr'

GR_for_pcls_savepath = 'GR_for_pcls'

sGR_for_pcls_savepath = 'sGR_for_pcls'

moa_gmt_savepath = 'moas.gmt'

corr_for_pcls_savepath = 'sGR_for_pcls_pearson_corr'

corr_rank_for_pcls_savepath = 'sGR_for_pcls_pearson_corr_rank'

% loocv outputs

corr_for_loo_cmpd_savepath = 'sGR_for_loo_cmpd_pearson_corr'

corr_for_loo_cmpd_to_remaining_reference_set_savepath = 'sGR_for_loo_cmpd_to_remaining_reference_set_pearson_corr'

col_meta_loo_cmpd_savepath = 'col_meta_loo_cmpd.txt'

% ## Load reference set MOA annotation sheet

annot = xls2table(reference_set_annot_path,1,true);

size(annot)

annot.pert_id

% ## Load sGR gct

if load_and_save_GR
    gr = parse_gctx(fullfile(datadir,'GR_filename'));
end

sgr = parse_gctx(fullfile(datadir,sGR_filename));

% ## Generate col_meta

[~,col_meta] = gct2meta(sgr);

headt(col_meta)

% ## Target description field is annotated compound MOA and is what will be used to describe each PCL cluster

annot.pcl_desc = annot.target_description;

fields = {'pert_id','pcl_desc'};

pcl_annot = [annot(:,fields)];

headt(col_meta)

size(col_meta)
[col_meta, ileft, ] = outerjoin(col_meta,pcl_annot,'Keys','pert_id','RightVariables','pcl_desc','Type','left','MergeKeys',true);

col_meta.target_description = col_meta.pcl_desc;
col_meta.proj_broad_id = strcat(col_meta.project_id, ':', col_meta.broad_id);

% re-sort col_meta data by original indices, caution: join operations in Matlab re-sort tables by their key values
[sorted_ileft, ileft_index] = sort(ileft);
col_meta = col_meta(ileft_index,:);
size(col_meta)

headt(col_meta)

% Save column metadata for all reference set and experimental conditions

wtable(col_meta,fullfile(wkdir, col_meta_savepath))

%col_meta = sortrows(col_meta,{'pcl_desc','pert_id','pert_dose'});

% Annotate sGR GCT with col_meta target_description (MOA) if absent or getting updated according to reference set annotation file

sgr = annotate_ds(sgr, table2struct(col_meta(:,{'cid','proj_broad_id','target_description'})),'dim','column','keyfield','cid');

% ## Keep only entries that are annotated with MOA (reference set CGI profiles)

tabulate(sort(col_meta.pcl_desc))

idx = cellfun(@isempty, col_meta.pcl_desc);
sum(idx)
col_meta(idx,:) = [];

% save column metadata for all reference set conditions (one row per CGI profile)

wtable(col_meta, fullfile(wkdir, col_meta_reference_set_savepath))

% ## Unwrap annotated MOA for compounds with multiple MOAs separated using "|" separator into multiple rows

size(col_meta)
col_meta_for_pcls = struct2table(unwrap_table(table2struct(col_meta),'pcl_desc','|'));
size(col_meta_for_pcls)

% Save column metadata for all reference set conditions as input for PCL clustering (one row per CGI profile-annotated MOA
% multiple rows for the same CGI profile if has multiple annotated MOAs

wtable(col_meta_for_pcls, fullfile(wkdir, col_meta_reference_set_for_pcls_savepath))

col_meta_for_pcls.pcl_desc = any2str(col_meta_for_pcls.pcl_desc);

% ## Calculate Pearson correlation between all reference set and experimental compound CGI profiles using sGR

sgr_corr = ds_corr(sgr, 'type', 'pearson');

mkgctx(fullfile(wkdir, corr_savepath),sgr_corr)

if save_gct
    mkgct(fullfile(wkdir, corr_savepath),sgr_corr,'precision',4)
end

% ## Subset gcts to reference set CGI profiles only as input for PCL clustering 

length(sgr.cid)

%gr.cid(1:10)
sgr.cid(1:10)

if load_and_save_GR
    gr_for_pcls = ds_slice(gr,'cid',unique(col_meta_for_pcls.cid,'stable'));
end

sgr_for_pcls = ds_slice(sgr,'cid',unique(col_meta_for_pcls.cid,'stable'));

if load_and_save_GR
    mkgctx(fullfile(wkdir, GR_for_pcls_savepath),gr)
end
mkgctx(fullfile(wkdir, sGR_for_pcls_savepath),sgr_for_pcls)

length(sgr_for_pcls.cid)

sgr_for_pcls.cid(1:10)

% ## Create gmt file with CGI profile - MOA membership

moas = tbl2gmt(table2struct(sortrows(col_meta_for_pcls, {'pcl_desc'})),'group_field','pcl_desc','desc_field','pcl_desc','member_field','cid')

% ## Remove any ambiguous MOA annotations

idx = ismember({moas.head}, filter_moa_annotations);
sum(idx)
moas(idx) = []

mkgmt(fullfile(wkdir, moa_gmt_savepath), moas)

% ## Calculate Pearson correlation and average rank of correlation (mutual nearest-neighbors) between all reference set CGI profiles using sGR 

sgr_corr_for_pcls = ds_corr(sgr_for_pcls, 'type', 'pearson');

sgr_corr_rank_for_pcls = sgr_corr_for_pcls;

% for each reference set condition, rank all reference set conditions from highest similarity/Pearson correlation to lowest similarity/Pearson correlation
mat = rankorder(sgr_corr_for_pcls.mat,'dim','row','direc','descend');

% in order to symmetrize matrix, calculate as the average rank between each of reference set conditions
sgr_corr_rank_for_pcls.mat = (mat+mat')/2;

mkgctx(fullfile(wkdir, corr_for_pcls_savepath),sgr_corr_for_pcls)
mkgctx(fullfile(wkdir, corr_rank_for_pcls_savepath),sgr_corr_rank_for_pcls)

if save_gct
    mkgct(fullfile(wkdir, corr_for_pcls_savepath),sgr_corr_for_pcls,'precision',4)
    mkgct(fullfile(wkdir, corr_rank_for_pcls_savepath),sgr_corr_rank_for_pcls,'precision',4)
end

% ## Save list of reference set compounds (pert_ids) for referencing during LOOCV

unique_reference_set_cmpds_list = unique(col_meta.pert_id);

length(unique_reference_set_cmpds_list)

unique_reference_set_cmpds_idx = (1:length(unique_reference_set_cmpds_list))';

unique_reference_set_cmpds_tbl = table(unique_reference_set_cmpds_idx, unique_reference_set_cmpds_list, 'VariableNames', {'reference_set_cmpd_idx', 'reference_set_cmpd'}); 

headt(unique_reference_set_cmpds_tbl)

wtable(unique_reference_set_cmpds_tbl, fullfile(wkdir,'reference_set_pert_ids_tbl_for_loocv.txt'))

% ## For each reference set Pert ID (compound) leave all doses of it out(Broad ID/Proj-Broad ID/CID) from the input correlation matrices and correlation rank matrices for PCL cluster construction 

if prepare_loocv

    all_reference_set_cids = unique(col_meta_for_pcls.cid, 'stable');

    length(all_reference_set_cids)
    
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

        col_meta_for_loo_cmpd = col_meta_for_pcls(ismember(col_meta_for_pcls.pert_id, unique_reference_set_cmpds_list(i)), :);

        assert(size(col_meta_for_loo_cmpd, 1) > 1)

        col_meta_for_pcls_remaining_reference_set = col_meta_for_pcls(~ismember(col_meta_for_pcls.pert_id, unique_reference_set_cmpds_list(i)), :);

        assert(size(col_meta_for_pcls_remaining_reference_set, 1) < size(col_meta_for_pcls, 1))

        wtable(col_meta_for_pcls_remaining_reference_set, fullfile(loo_wkdir, col_meta_reference_set_for_pcls_savepath));

        wtable(col_meta_for_loo_cmpd, fullfile(loo_wkdir, col_meta_loo_cmpd_savepath));




        remaining_reference_set_cids = unique(col_meta_for_pcls_remaining_reference_set.cid,'stable');

        loo_cmpd_cids = unique(col_meta_for_loo_cmpd.cid,'stable');

        assert(length(remaining_reference_set_cids) < length(all_reference_set_cids))

        assert(length(loo_cmpd_cids) > 0)

        sgr_corr_remaining_reference_set = ds_slice(sgr_corr_for_pcls, 'cid',remaining_reference_set_cids, 'rid', remaining_reference_set_cids);

        sgr_corr_loo_cmpd = ds_slice(sgr_corr_for_pcls, 'cid',loo_cmpd_cids, 'rid', loo_cmpd_cids);

        sgr_corr_loo_cmpd_to_remaining_reference_set = ds_slice(sgr_corr_for_pcls, 'rid',remaining_reference_set_cids, 'cid', loo_cmpd_cids);

        assert(size(sgr_corr_loo_cmpd_to_remaining_reference_set.mat, 1) < size(sgr_corr_for_pcls.mat, 1))

        assert(size(sgr_corr_loo_cmpd.mat, 1) > 0)


        sgr_corr_rank_remaining_reference_set = sgr_corr_remaining_reference_set;
        mat = rankorder(sgr_corr_remaining_reference_set.mat,'dim','row','direc','descend');
        sgr_corr_rank_remaining_reference_set.mat = (mat+mat')/2;

        mkgctx(fullfile(loo_wkdir, corr_for_pcls_savepath),sgr_corr_remaining_reference_set);
        mkgctx(fullfile(loo_wkdir, corr_rank_for_pcls_savepath),sgr_corr_rank_remaining_reference_set);

        mkgctx(fullfile(loo_wkdir, corr_for_loo_cmpd_savepath),sgr_corr_loo_cmpd);
        mkgctx(fullfile(loo_wkdir, corr_for_loo_cmpd_to_remaining_reference_set_savepath),sgr_corr_loo_cmpd_to_remaining_reference_set);

        moas_remaining_reference_set = tbl2gmt(table2struct(sortrows(col_meta_for_pcls_remaining_reference_set, {'pcl_desc'})),'group_field','pcl_desc','desc_field','pcl_desc','member_field','cid');

        idx = ismember({moas_remaining_reference_set.head}, filter_moa_annotations);
        sum(idx);
        moas_remaining_reference_set(idx) = [];

        mkgmt(fullfile(loo_wkdir, moa_gmt_savepath), moas_remaining_reference_set);
    end
end

