% # Start

datadir = '../data';
wkdir = '../results';

mk_cd_dir(wkdir, true);
%imatlab_export_fig('print-png')

% inputs

kabx_annot_path = '../data/kabx_annotations.xlsx'

GR_filename = 'GR_kabx_gsk_brd4310_n10819x340.gctx'

sGR_filename = 'sGR_kabx_gsk_brd4310_n10819x340.gctx'

filter_moa_annotations = {'NaN','unknown','whole cell only'}

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

% ## Load KABX MOA annotation sheet

annot = xls2table(kabx_annot_path,1,true);

size(annot)

annot.pert_id

% ## Load sGR gct

%gr = parse_gctx(fullfile(datadir,'GR_filename'));

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

%mkgct(fullfile(wkdir, corr_savepath),grzs_corr,'precision',4)

% ## Subset gcts to KABX dsCGI profiles only as input for PCL clustering 

length(grzs.cid)

%gr.cid(1:10)
grzs.cid(1:10)

%gr_for_pcls = ds_slice(gr,'cid',unique(col_meta_for_pcls.cid,'stable'));
grzs_for_pcls = ds_slice(grzs,'cid',unique(col_meta_for_pcls.cid,'stable'));

%mkgctx(fullfile(wkdir, GR_for_pcls_savepath),gr)
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

%mkgct(fullfile(wkdir, corr_for_pcls_savepath),grzs_corr_for_pcls,'precision',4)
%mkgct(fullfile(wkdir, corr_rank_for_pcls_savepath),grzs_corr_rank_for_pcls,'precision',4)