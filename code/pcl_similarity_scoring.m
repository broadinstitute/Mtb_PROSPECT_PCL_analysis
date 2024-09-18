function run_pcls_validation_no_loo_scoring_exclude_self(prefix,project_id,pcls_moa_path,c_path,c_rank_path,col_meta_path,cp_annot_path,fgr_path,wkdir)

% RUN_PCLS_VALIDATION(PREFIX,PROJECT_ID,PCLS_MOA,DS_CORR,DS_CORR_RANK,COL_META,CP_ANNOT,FGR,WKDIR)
% Set initial variables
% prefix = 'known_moas_all';
% project_id = 'kabx2';
% 
% datadir = '/idi/cgtb/morzech/idmp/screen4wk_screen5wk_kabx2_tbda1/analysis/pcls';
% gmtdir = fullfile(datadir, prefix, 'iteration_final');
% wkdir = fullfile(datadir, prefix, 'iteration_final',['pcls_for_',project_id,'_final4']);
mk_cd_dir(wkdir, true)

% Load final pcls
if ischar(pcls_moa_path)
	pcls_moa = parse_gmt(pcls_moa_path);
end
% Remove groups with fewer than 2 profiles
num_pcls_ini = numel(pcls_moa);
disp(sprintf('The initial number of PCLs is %d', num_pcls_ini))
pcls_moa([pcls_moa.len]<2) = [];
num_pcls_final = numel(pcls_moa);
disp(sprintf('The remaining number of PCLs is %d', num_pcls_final))

% Load compound annotation
if isstr(cp_annot_path)
    [d,f,e] = fileparts(cp_annot_path);
    if ismember(e,{'.txt'})==1
		cp_annot = rtable(cp_annot_path);
    elseif ismember(e,{'.xls','.xlsx'})==1
		cp_annot = xls2table(cp_annot_path,1,true);
    else
		error('Unknown file format for compound annotation! Should be txt, xls, or xlsx')
    end
end

% Load table with fractions of GR
if isstr(fgr_path)
    fgr = rtable(fgr_path);
end

if ismember('double',unique(class(fgr.x_median_total_count)))~=1
    fgr.x_median_total_count = cellfun(@str2double, fgr.x_median_total_count);
end

if ismember('double',unique(class(fgr.frac_gr_le0_30)))~=1
    fgr.frac_gr_le0_30 = cellfun(@str2double, fgr.frac_gr_le0_30);
end

% Load column metadata
if isstr(col_meta_path)
	col_meta = rtable(col_meta_path);
end
if ismember('double',unique(class(col_meta.x_median_total_count)))~=1
    col_meta.x_median_total_count = cellfun(@str2double, col_meta.x_median_total_count);
end

% Add sum_gr_le0_30 and frac_gr_le0_30 to col_meta if not already there
if ismember('sum_gr_le0_30', col_meta.Properties.VariableNames)==0
	col_meta = join(col_meta, fgr,'Keys','cid','RightVariables','sum_gr_le0_30');
end

if ismember('frac_gr_le0_30', col_meta.Properties.VariableNames)==0
	col_meta = join(col_meta, fgr,'Keys','cid','RightVariables','frac_gr_le0_30');
end

% Extract entries with a given project_id and not well-killers
if ~isempty(project_id)
	rids = col_meta.cid(ismember(col_meta.project_id, project_id) & col_meta.frac_gr_le0_30 < 0.65);
else
	%rids = col_meta.cid(col_meta.frac_gr_le0_30 < 0.65);
    rids = col_meta.cid;
end

% Load Pearson correlation for the combined dataset
if isstr(c_path)
	c = parse_gctx(c_path);
end
c = annotate_ds(c,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','column','keyfield','cid');
c = annotate_ds(c,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','row','keyfield','cid');

% Load correlation-based rank
if isstr(c_rank_path)
	c_rank = parse_gctx(c_rank_path);
end
c_rank = annotate_ds(c_rank,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','column','keyfield','cid');
c_rank = annotate_ds(c_rank,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','row','keyfield','cid');

% Make sure that c and c_rank are sorted in the same way
c_rank = ds_slice(c_rank,'cid',c.cid,'rid',c.rid);

assert(all(diag(c.mat) == 1))

assert(all(diag(c_rank.mat) == 1))

% Set the diagonal of the matrix to all NaN values
c.mat(eye(size(c.mat))==1) = NaN;

c_rank.mat(eye(size(c_rank.mat))==1) = NaN;

assert(all(isnan(diag(c.mat))))

assert(all(isnan(diag(c_rank.mat))))

% Add pert_id
c.rdesc(:, c.rdict('pert_id')) = c.rdesc(:, c.rdict('broad_id'));
idx = find_match(c.rdesc(:, c.rdict('broad_id')),'BRD-',true);
c.rdesc(idx, c.rdict('pert_id')) = cellfun(@(s) s(1:13), c.rdesc(idx, c.rdict('broad_id')), 'uni', 0);

%% Check if all reference compounds are in the correlation matrix
setdiff(cat(1,pcls_moa.entry),c.cid)

% Slice c and c_rank to contain all profiles in rows and profiles from PCLs in columns
ridx = ismember(c.rid, rids);
sum(ridx)
cidx = ismember(c.cid, cat(1,pcls_moa.entry));
sum(cidx)

c = ds_slice(c,'ridx',ridx,'cidx',cidx);
c_rank = ds_slice(c_rank,'ridx',ridx,'cidx',cidx);

% Save gctx files
disp('Saving gctx files')
mkgctx(fullfile(wkdir,'ds_corr.gctx'), c)
mkgctx(fullfile(wkdir,'ds_corr_rank.gctx'), c_rank)

% ## Run MoA discovery using MoA-based PCLs with overlap
assert(isequal(c.cid, c_rank.cid),'Columns are not sorted identically')
assert(isequal(c.rid, c_rank.rid),'Rows are not sorted identically')

% pre-sort indices by broad_id
disp('Sorting indices by broad_id')
issorted(c.rdesc(:, c.rdict('broad_id')))
[broad_id_sort,broad_id_sort_idx] = sort(c.rdesc(:, c.rdict('broad_id')));
issorted(c.rdesc(broad_id_sort_idx, c.rdict('broad_id')))
numel(broad_id_sort_idx)

% Initialize output structure
iteration = [1:numel(c.rid)]';
tbl = repmat({'empty'},numel(c.rid),1);
out = table2struct(table(iteration,tbl));
pcls = pcls_moa;
num_pcls = numel(pcls)

out_tbl = {}

% commented out to avoid parallel pool initialization issue
%p = parpool(feature('numcores'),'IdleTimeout',120);

% Added by AB on 06/12/2024
% Ensure parallel pool is cleaned up at the end of the script or on error
cleanupObj = onCleanup(@() delete(gcp('nocreate')));

%parfor jj = 1:numel(pcls)
for jj = 1:numel(pcls)

    loop_progress(jj,numel(pcls),100)

    tmp_tbl = {};
    
    tmp_tbl.rid = c.rid;
    tmp_tbl.broad_id = c.rdesc(:, c.rdict('broad_id'));
    tmp_tbl.pert_id = c.rdesc(:, c.rdict('pert_id'));
    tmp_tbl.pert_idose = c.rdesc(:, c.rdict('pert_idose'));
    tmp_tbl.pert_dose = c.rdesc(:, c.rdict('pert_dose'));
    tmp_tbl.x_median_total_count = c.rdesc(:, c.rdict('x_median_total_count'));

    num_rids = numel(c.rid);

    %pcls(jj).head

    tmp_tbl.pcl_id = repmat({pcls(jj).head}, num_rids, 1);

    tmp_tbl.pcl_desc = repmat({pcls(jj).desc}, num_rids, 1);

    tmp_tbl.pcl_size = repmat(pcls(jj).len, num_rids, 1);

    cidx = ismember(c.cid,pcls(jj).entry);

    tmp_tbl.pcl_size_actual = repmat(sum(cidx), num_rids, 1);

    tmp_tbl.pcl_broad_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('broad_id'))))}, num_rids, 1);
    tmp_tbl.pcl_pert_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('pert_id'))))}, num_rids, 1);
    tmp_tbl.num_pert_id_unique = repmat(numel(unique(c.cdesc(cidx, c.cdict('pert_id')))), num_rids, 1);
    tmp_tbl.median_corr = median(c.mat(:, cidx), 2, 'omitnan');
    tmp_tbl.median_corr_rank = median(c_rank.mat(:,cidx), 2, 'omitnan');

    tmp_tbl = struct2table(tmp_tbl);

    %tmp_tbl = sortrows(tmp_tbl,{'broad_id','median_corr_rank'},{'ascend','ascend'});
    
    out_tbl = [out_tbl;tmp_tbl];
    
end

size(out_tbl)

disp('All PCLs successfully looped through')

out_tbl.median_corr_rank_div_group_size = out_tbl.median_corr_rank./out_tbl.pcl_size_actual;

headt(out_tbl)

% Make gctx
out_tbl.rid_idx = grp2idx(out_tbl.rid);
out_tbl.pcl_id_idx = grp2idx(out_tbl.pcl_id);

% Find indices
[~,ridx] = unique(out_tbl.pcl_id_idx);
rid = out_tbl.pcl_id(ridx);
[~,cidx] = unique(out_tbl.rid_idx);
cid = out_tbl.rid(cidx);

% Create metadata
col_meta = out_tbl(cidx,{'rid','broad_id','pert_id','pert_idose','pert_dose'});
row_meta = out_tbl(ridx,{'pcl_id','pcl_desc','pcl_size','pcl_broad_id','pcl_pert_id','num_pert_id_unique'});

% Create an empty matrix
mat = nan(numel(rid),numel(cid));
ind = sub2ind(size(mat),out_tbl.pcl_id_idx,out_tbl.rid_idx);

ds = mkgctstruct(mat,'rid',rid','cid',cid);
ds = annotate_ds(ds,table2struct(col_meta),'dim','column','keyfield','rid');
ds = annotate_ds(ds,table2struct(row_meta),'dim','row','keyfield','pcl_id');

% Create gcts with various data
fields = {'median_corr','median_corr_rank','median_corr_rank_div_group_size'};

for ii = 1:numel(fields)
    disp(fields{ii})
    ds_tmp = ds;
    ds_tmp.mat(ind) = out_tbl.(fields{ii});

    % Verify
    rid1 = ds_tmp.rid(2);
    cid1 = ds_tmp.cid(1);
    corr1 = ds_tmp.mat(2,1);
    idx = ismember(out_tbl.rid,cid1)&ismember(out_tbl.pcl_id,rid1);
    assert(sum(idx)>0,'Rows or columns of the gct structure do not match');
    assert(out_tbl.(fields{ii})(idx)==corr1,'Value in the matrix do not match the one in the table');
    mkgctx(fullfile(wkdir,['pcl_',fields{ii},'.gctx']), ds_tmp)
end

% Save the output table
%if ~isempty(project_id)
%    wtable(out_tbl,fullfile(wkdir,[project_id, '_correlation_to_selected_pcls.txt']));
%else
%    wtable(out_tbl,fullfile(wkdir,'correlation_to_selected_pcls.txt'));
%end

