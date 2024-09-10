% # Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, true);
%imatlab_export_fig('print-png')

% loocv inputs

prepare_loocv = true

demo_loocv = true

demo_loocv_number_or_list = 'number' % 'number' or 'list'

demo_loocv_number_cmpds = 2

demo_loocv_list_cmpds = {'BRD-K04804440','BRD-K01507359','BRD-K87202646','BRD-K59853741', 'BRD-K27302037'} % Ciprofloxacin, Rifampin, Isoniazid, Q203, Thioacetazone

results_subdir_prefix = 'loocv_pcls/leave_out_cmpd_'

% inputs

kabx_annot_path = '../data/kabx_annotations.xlsx'

GR_filename = 'GR_kabx_gsk_brd4310_n10819x340.gctx'

sGR_filename = 'sGR_kabx_gsk_brd4310_n10819x340.gctx'

filter_moa_annotations = {'NaN','unknown','whole cell only'}

load_and_save_GR = false

save_gct = false % greater memory consumption than .gctx but can be opened in Morpheus, https://software.broadinstitute.org/morpheus

% outputs

col_meta_savepath = 'col_meta.txt'

col_meta_kabx_savepath = 'col_meta_kabx.txt'

col_meta_kabx_for_pcls_savepath = 'col_meta_kabx_for_pcls.txt'

corr_savepath = 'sGR_kabx_gsk_brd4310_pearson_corr'

GR_for_pcls_savepath = 'GR_for_pcls'

sGR_for_pcls_savepath = 'sGR_for_pcls'

moa_gmt_savepath = 'moas.gmt'

corr_for_pcls_savepath = 'sGR_for_pcls_pearson_corr'

corr_rank_for_pcls_savepath = 'sGR_for_pcls_pearson_corr_rank'

% loocv outputs

corr_for_loo_cmpd_savepath = 'sGR_for_loo_cmpd_pearson_corr'

corr_for_loo_cmpd_to_remaining_kabx_savepath = 'sGR_for_loo_cmpd_to_remaining_kabx_pearson_corr'

col_meta_loo_cmpd_savepath = 'col_meta_loo_cmpd.txt'

% ## Load KABX MOA annotation sheet

annot = xls2table(kabx_annot_path,1,true);

size(annot)

annot.pert_id

% ## Load sGR gct

if load_and_save_GR
    gr = parse_gctx(fullfile(datadir,'GR_filename'));
end

grzs = parse_gctx(fullfile(datadir,sGR_filename));

% ## Generate col_meta

[~,col_meta] = gct2meta(grzs);

headt(col_meta)

% ## Target description field is annotated compound MOA and is what will be used to describe each PCL

annot.pcl_desc = annot.target_description;

fields = {'pert_id','pcl_desc'};

pcl_annot = [annot(:,fields)];

headt(col_meta)

size(col_meta)
[col_meta, ileft, ] = outerjoin(col_meta,pcl_annot,'Keys','pert_id','RightVariables','pcl_desc','Type','left','MergeKeys',true);

% re-sort col_meta data by original indices, caution: join operations in Matlab re-sort tables by their key values
[sorted_ileft, ileft_index] = sort(ileft);
col_meta = col_meta(ileft_index,:);
size(col_meta)

headt(col_meta)

% Save column metadata for all KABX and experimental treatments

wtable(col_meta,fullfile(wkdir, col_meta_savepath))

%col_meta = sortrows(col_meta,{'pcl_desc','pert_id','pert_dose'});

% ## Keep only entries that are annotated with MOA (KABX dsCGI profiles)

tabulate(sort(col_meta.pcl_desc))

idx = cellfun(@isempty, col_meta.pcl_desc);
sum(idx)
col_meta(idx,:) = [];

% save column metadata for all KABX treatments (one row per dsCGI profile)

wtable(col_meta, fullfile(wkdir, col_meta_kabx_savepath))

% ## Unwrap annotated MOA for compounds with multiple MOAs separated using "|" separator into multiple rows

size(col_meta)
col_meta_for_pcls = struct2table(unwrap_table(table2struct(col_meta),'pcl_desc','|'));
size(col_meta_for_pcls)

% Save column metadata for all KABX treatments as input for PCL clustering (one row per dsCGI profile-annotated MOA
% multiple rows for the same dsCGI profile if has multiple annotated MOAs

wtable(col_meta_for_pcls, fullfile(wkdir, col_meta_kabx_for_pcls_savepath))

col_meta_for_pcls.pcl_desc = any2str(col_meta_for_pcls.pcl_desc);

% ## Calculate Pearson correlation between all KABX and experimental compound dsCGI profiles using sGR

grzs_corr = ds_corr(grzs, 'type', 'pearson');

mkgctx(fullfile(wkdir, corr_savepath),grzs_corr)

if save_gct
    mkgct(fullfile(wkdir, corr_savepath),grzs_corr,'precision',4)
end

% ## Subset gcts to KABX dsCGI profiles only as input for PCL clustering 

length(grzs.cid)

%gr.cid(1:10)
grzs.cid(1:10)

if load_and_save_GR
    gr_for_pcls = ds_slice(gr,'cid',unique(col_meta_for_pcls.cid,'stable'));
end

grzs_for_pcls = ds_slice(grzs,'cid',unique(col_meta_for_pcls.cid,'stable'));

if load_and_save_GR
    mkgctx(fullfile(wkdir, GR_for_pcls_savepath),gr)
end
mkgctx(fullfile(wkdir, sGR_for_pcls_savepath),grzs_for_pcls)

length(grzs_for_pcls.cid)

grzs_for_pcls.cid(1:10)

% ## Create gmt file with dsCGI profile - MOA membership

moas = tbl2gmt(table2struct(sortrows(col_meta_for_pcls, {'pcl_desc'})),'group_field','pcl_desc','desc_field','pcl_desc','member_field','cid')

% ## Remove any ambiguous MOA annotations

idx = ismember({moas.head}, filter_moa_annotations);
sum(idx)
moas(idx) = []

mkgmt(fullfile(wkdir, moa_gmt_savepath), moas)

% ## Calculate Pearson correlation and average rank of correlation (mutual nearest-neighbors) between all KABX dsCGI profiles using sGR 

grzs_corr_for_pcls = ds_corr(grzs_for_pcls, 'type', 'pearson');

grzs_corr_rank_for_pcls = grzs_corr_for_pcls;

% for each KABX treatment, rank all KABX treatments from highest similarity/Pearson correlation to lowest similarity/Pearson correlation
mat = rankorder(grzs_corr_for_pcls.mat,'dim','row','direc','descend');

% in order to symmetrize matrix, calculate as the average rank between each of KABX treatments
grzs_corr_rank_for_pcls.mat = (mat+mat')/2;

mkgctx(fullfile(wkdir, corr_for_pcls_savepath),grzs_corr_for_pcls)
mkgctx(fullfile(wkdir, corr_rank_for_pcls_savepath),grzs_corr_rank_for_pcls)

if save_gct
    mkgct(fullfile(wkdir, corr_for_pcls_savepath),grzs_corr_for_pcls,'precision',4)
    mkgct(fullfile(wkdir, corr_rank_for_pcls_savepath),grzs_corr_rank_for_pcls,'precision',4)
end

% ## Save list of KABX compounds (pert_ids) for referencing during LOOCV

unique_kabx_cmpds_list = unique(col_meta.pert_id);

length(unique_kabx_cmpds_list)

unique_kabx_cmpds_idx = (1:length(unique_kabx_cmpds_list))';

unique_kabx_cmpds_tbl = table(unique_kabx_cmpds_idx, unique_kabx_cmpds_list, 'VariableNames', {'kabx_cmpd_idx', 'kabx_cmpd'}); 

headt(unique_kabx_cmpds_tbl)

wtable(unique_kabx_cmpds_tbl, fullfile(wkdir,'kabx_pert_ids_tbl_for_loocv.txt'))

% ## For each KABX Pert ID (compound) leave all doses of it out(Broad ID/Proj-Broad ID/CID) from the input correlation matrices and correlation rank matrices for PCL construction 

if prepare_loocv

    all_kabx_cids = unique(col_meta_for_pcls.cid, 'stable');

    length(all_kabx_cids)
    
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

        col_meta_for_loo_cmpd = col_meta_for_pcls(ismember(col_meta_for_pcls.pert_id, unique_kabx_cmpds_list(i)), :);

        assert(size(col_meta_for_loo_cmpd, 1) > 1)

        col_meta_for_pcls_remaining_kabx = col_meta_for_pcls(~ismember(col_meta_for_pcls.pert_id, unique_kabx_cmpds_list(i)), :);

        assert(size(col_meta_for_pcls_remaining_kabx, 1) < size(col_meta_for_pcls, 1))

        wtable(col_meta_for_pcls_remaining_kabx, fullfile(loo_wkdir, col_meta_kabx_for_pcls_savepath));

        wtable(col_meta_for_loo_cmpd, fullfile(loo_wkdir, col_meta_loo_cmpd_savepath));




        remaining_kabx_cids = unique(col_meta_for_pcls_remaining_kabx.cid,'stable');

        loo_cmpd_cids = unique(col_meta_for_loo_cmpd.cid,'stable');

        assert(length(remaining_kabx_cids) < length(all_kabx_cids))

        assert(length(loo_cmpd_cids) > 0)

        grzs_corr_remaining_kabx = ds_slice(grzs_corr_for_pcls, 'cid',remaining_kabx_cids, 'rid', remaining_kabx_cids);

        grzs_corr_loo_cmpd = ds_slice(grzs_corr_for_pcls, 'cid',loo_cmpd_cids, 'rid', loo_cmpd_cids);

        grzs_corr_loo_cmpd_to_remaining_kabx = ds_slice(grzs_corr_for_pcls, 'rid',remaining_kabx_cids, 'cid', loo_cmpd_cids);

        assert(size(grzs_corr_loo_cmpd_to_remaining_kabx.mat, 1) < size(grzs_corr_for_pcls.mat, 1))

        assert(size(grzs_corr_loo_cmpd.mat, 1) > 0)


        grzs_corr_rank_remaining_kabx = grzs_corr_remaining_kabx;
        mat = rankorder(grzs_corr_remaining_kabx.mat,'dim','row','direc','descend');
        grzs_corr_rank_remaining_kabx.mat = (mat+mat')/2;

        mkgctx(fullfile(loo_wkdir, corr_for_pcls_savepath),grzs_corr_remaining_kabx);
        mkgctx(fullfile(loo_wkdir, corr_rank_for_pcls_savepath),grzs_corr_rank_remaining_kabx);

        mkgctx(fullfile(loo_wkdir, corr_for_loo_cmpd_savepath),grzs_corr_loo_cmpd);
        mkgctx(fullfile(loo_wkdir, corr_for_loo_cmpd_to_remaining_kabx_savepath),grzs_corr_loo_cmpd_to_remaining_kabx);

        moas_remaining_kabx = tbl2gmt(table2struct(sortrows(col_meta_for_pcls_remaining_kabx, {'pcl_desc'})),'group_field','pcl_desc','desc_field','pcl_desc','member_field','cid');

        idx = ismember({moas_remaining_kabx.head}, filter_moa_annotations);
        sum(idx);
        moas_remaining_kabx(idx) = [];

        mkgmt(fullfile(loo_wkdir, moa_gmt_savepath), moas_remaining_kabx);
    end
end

