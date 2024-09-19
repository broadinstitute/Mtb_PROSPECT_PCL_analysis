function pcl_similarity_scoring(clusters_path,c_path,c_rank_path,col_meta_path,fgr_path,wkdir,min_clust_size)

% PCL_SIMILARITY_SCORING(CLUSTERS_PATH,DS_CORR,DS_CORR_RANK,COL_META,FGR,WKDIR,MIN_CLUST_SIZE)
% Set initial variables
% 
% datadir = '/idi/cgtb/morzech/idmp/screen4wk_screen5wk_kabx2_tbda1/analysis/pcls';
% gmtdir = fullfile(datadir, prefix, 'iteration_final');
% wkdir = fullfile(datadir, prefix, 'iteration_final',['pcls_for_',project_id,'_final4']);
mk_cd_dir(wkdir, false)

% Load clusters from spectral clustering
if ischar(clusters_path)
	clusters = parse_gmt(clusters_path);
end

% Set minimum cluster size if not provided in the input
if isempty(min_clust_size)
	min_clust_size = 2; 
end

% Remove groups with fewer than min_clust_size profiles
num_clusters_ini = numel(clusters);
disp(sprintf('The initial number of clusters is %d', num_clusters_ini))
disp(sprintf('Removing clusters with size less than %d', min_clust_size))
clusters([clusters.len]<min_clust_size) = [];
num_clusters_final = numel(clusters);
disp(sprintf('The remaining number of clusters is %d', num_clusters_final))

% Load compound annotation
% was not used for anything
%if isstr(cp_annot_path)
%    [d,f,e] = fileparts(cp_annot_path);
%    if ismember(e,{'.txt'})==1
%		cp_annot = rtable(cp_annot_path);
%    elseif ismember(e,{'.xls','.xlsx'})==1
%		cp_annot = xls2table(cp_annot_path,1,true);
%    else
%		error('Unknown file format for compound annotation! Should be txt, xls, or xlsx')
%    end
%end

%% Load table with fractions of GR
%if isstr(fgr_path)
%    fgr = rtable(fgr_path);
%end

%if ismember('double',unique(class(fgr.x_median_total_count)))~=1
%    fgr.x_median_total_count = cellfun(@str2double, fgr.x_median_total_count);
%end

%if ismember('double',unique(class(fgr.frac_gr_le0_30)))~=1
%    fgr.frac_gr_le0_30 = cellfun(@str2double, fgr.frac_gr_le0_30);
%end

% Load column metadata
if isstr(col_meta_path)
	col_meta = rtable(col_meta_path);
end

%if ismember('double',unique(class(col_meta.x_median_total_count)))~=1
%    col_meta.x_median_total_count = cellfun(@str2double, col_meta.x_median_total_count);
%end

%% Add sum_gr_le0_30 and frac_gr_le0_30 to col_meta if not already there
%if ismember('sum_gr_le0_30', col_meta.Properties.VariableNames)==0
%	col_meta = join(col_meta, fgr,'Keys','cid','RightVariables','sum_gr_le0_30');
%end

%if ismember('frac_gr_le0_30', col_meta.Properties.VariableNames)==0
%	col_meta = join(col_meta, fgr,'Keys','cid','RightVariables','frac_gr_le0_30');
%end

% Extract entries with a given project_id and not well-killers
% have commented out/disabled functionality for using fraction of strains with GR below a threshold to define and exclude well-killers
%if ~isempty(project_id)
%	rids = col_meta.cid(ismember(col_meta.project_id, project_id));
	%rids = col_meta.cid(ismember(col_meta.project_id, project_id) & col_meta.frac_gr_le0_30 < 0.65);
%else
	%rids = col_meta.cid(col_meta.frac_gr_le0_30 < 0.65);
%    rids = col_meta.cid;
%end

rids = col_meta.cid;

% Load Pearson correlation for the combined dataset
if isstr(c_path)
	c = parse_gctx(c_path);
end
%c = annotate_ds(c,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','column','keyfield','cid');
%c = annotate_ds(c,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','row','keyfield','cid');

% Load correlation-based rank
if isstr(c_rank_path)
	c_rank = parse_gctx(c_rank_path);
end
%c_rank = annotate_ds(c_rank,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','column','keyfield','cid');
%c_rank = annotate_ds(c_rank,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','row','keyfield','cid');

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
% Already includes pert_id
%c.rdesc(:, c.rdict('pert_id')) = c.rdesc(:, c.rdict('broad_id'));
%idx = find_match(c.rdesc(:, c.rdict('broad_id')),'BRD-',true);
%c.rdesc(idx, c.rdict('pert_id')) = cellfun(@(s) s(1:13), c.rdesc(idx, c.rdict('broad_id')), 'uni', 0);

%% Check if all reference compounds are in the correlation matrix
setdiff(cat(1,clusters.entry),c.cid)

% Slice c and c_rank to contain all profiles in rows and profiles from clusters in columns
ridx = ismember(c.rid, rids);
sum(ridx)
cidx = ismember(c.cid, cat(1,clusters.entry));
sum(cidx)

c = ds_slice(c,'ridx',ridx,'cidx',cidx);
c_rank = ds_slice(c_rank,'ridx',ridx,'cidx',cidx);

% Save gctx files
disp('Saving gctx files')
mkgctx(fullfile(wkdir,'ds_corr.gctx'), c)
mkgctx(fullfile(wkdir,'ds_corr_rank.gctx'), c_rank)

% ## Run MoA discovery using MoA-based clusters with overlap
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
clusters = clusters;
num_clusters = numel(clusters)

out_tbl = {}

% Ensure parallel pool is cleaned up at the end of the script or on error
cleanupObj = onCleanup(@() delete(gcp('nocreate')));

% commented out to avoid parallel pool initialization issue
%p = parpool(feature('numcores'),'IdleTimeout',120);

%parfor jj = 1:numel(clusters)
for jj = 1:numel(clusters)

    loop_progress(jj,numel(clusters),100)

    tmp_tbl = {};
    
    tmp_tbl.rid = c.rid;
    tmp_tbl.proj_broad_id = c.rdesc(:, c.rdict('proj_broad_id'));
    tmp_tbl.broad_id = c.rdesc(:, c.rdict('broad_id'));
    tmp_tbl.pert_id = c.rdesc(:, c.rdict('pert_id'));
    tmp_tbl.pert_idose = c.rdesc(:, c.rdict('pert_idose_tickl'));
    tmp_tbl.pert_dose = c.rdesc(:, c.rdict('pert_dose'));
    %tmp_tbl.x_median_total_count = c.rdesc(:, c.rdict('x_median_total_count'));

    num_rids = numel(c.rid);

    %clusters(jj).head

    tmp_tbl.cluster_id = repmat({clusters(jj).head}, num_rids, 1);

    tmp_tbl.cluster_desc = repmat({clusters(jj).desc}, num_rids, 1);

    tmp_tbl.cluster_size = repmat(clusters(jj).len, num_rids, 1);

    cidx = ismember(c.cid,clusters(jj).entry);

    tmp_tbl.cluster_size_actual = repmat(sum(cidx), num_rids, 1);

    tmp_tbl.cluster_proj_broad_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('proj_broad_id'))))}, num_rids, 1);
    tmp_tbl.cluster_broad_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('broad_id'))))}, num_rids, 1);
    tmp_tbl.cluster_pert_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('pert_id'))))}, num_rids, 1);
    tmp_tbl.num_pert_id_unique = repmat(numel(unique(c.cdesc(cidx, c.cdict('pert_id')))), num_rids, 1);
    
    tmp_tbl.median_corr = median(c.mat(:, cidx), 2, 'omitnan');
    tmp_tbl.median_corr_rank = median(c_rank.mat(:,cidx), 2, 'omitnan');

    tmp_tbl = struct2table(tmp_tbl);

    %tmp_tbl = sortrows(tmp_tbl,{'broad_id','median_corr_rank'},{'ascend','ascend'});
    
    out_tbl = [out_tbl;tmp_tbl];
    
end

size(out_tbl)

disp('All clusters successfully looped through')

headt(out_tbl)

% Make gctx
out_tbl.rid_idx = grp2idx(out_tbl.rid);
out_tbl.cluster_id_idx = grp2idx(out_tbl.cluster_id);

% Find indices
[~,ridx] = unique(out_tbl.cluster_id_idx);
rid = out_tbl.cluster_id(ridx);
[~,cidx] = unique(out_tbl.rid_idx);
cid = out_tbl.rid(cidx);

% Create metadata

%col_meta = out_tbl(cidx,{'rid','proj_broad_id','broad_id','pert_id','pert_idose_tickl','pert_dose'});

row_meta = out_tbl(ridx,{'cluster_id','cluster_desc','cluster_size','cluster_proj_broad_id','cluster_broad_id','cluster_pert_id','num_pert_id_unique'});

% Create an empty matrix
mat = nan(numel(rid),numel(cid));
ind = sub2ind(size(mat),out_tbl.cluster_id_idx,out_tbl.rid_idx);

ds = mkgctstruct(mat,'rid',rid','cid',cid);
ds = annotate_ds(ds,table2struct(col_meta),'dim','column','keyfield','rid');
ds = annotate_ds(ds,table2struct(row_meta),'dim','row','keyfield','cluster_id');

% Create gctx files with various data
fields = {'median_corr','median_corr_rank'};

for ii = 1:numel(fields)
    disp(fields{ii})
    ds_tmp = ds;
    ds_tmp.mat(ind) = out_tbl.(fields{ii});

    % Verify
    rid1 = ds_tmp.rid(2);
    cid1 = ds_tmp.cid(1);
    corr1 = ds_tmp.mat(2,1);
    idx = ismember(out_tbl.rid,cid1)&ismember(out_tbl.cluster_id,rid1);
    assert(sum(idx)>0,'Rows or columns of the gct structure do not match');
    assert(out_tbl.(fields{ii})(idx)==corr1,'Value in the matrix do not match the one in the table');
    mkgctx(fullfile(wkdir,['clusters_',fields{ii},'.gctx']), ds_tmp)
end

% Section to check clusters for their most similar KABX dsCGI profile and define them as PCLs if the profiles are in-MOA (share same MOA as PCL)
% and exclude them as uninterpretable clusters otherwise

% Parse gctx with the cluster similarity score
g = glob(fullfile(wkdir,'cluster_median_corr_n*.gctx'))
s = parse_gctx(g{1});

% Select profiles with annotated MOA from KABX
cidx = ~cellfun(@isempty,s.cdesc(:,s.cdict('target_description')));
sum(cidx)

ss = ds_slice(s,'cidx',cidx)

col_meta_ss = cell2table([ss.cid,ss.cdesc],'VariableNames',['cid';ss.chd]);
col_meta_ss.target_description = any2str(col_meta_ss.target_description);
headt(col_meta_ss)
unique(col_meta_ss.target_description)

row_meta_ss = cell2table([ss.rid,ss.rdesc],'VariableNames',['rid';ss.rhd]);
row_meta_ss.pcl_desc = any2str(row_meta_ss.pcl_desc);
headt(row_meta_ss)
unique(row_meta_ss.pcl_desc)

[a,b] = ind2sub(size(ss.mat),1:numel(ss.mat));

ss_tbl = [col_meta_ss(b,:),row_meta_ss(a,:)];
size(ss_tbl)
ss_tbl.pcl_score = ss.mat(:);

ss_same_moa = ss;
ss_same_moa.mat = zeros(size(ss.mat));
target_desc = ss_same_moa.cdesc(:, ss_same_moa.cdict('target_description'));

for ii = 1:numel(ss_same_moa.rid)
    pcl_desc = ss_same_moa.rdesc(ii,ss_same_moa.rdict('pcl_desc'));
    %ss_same_moa.mat(ii,:) = ismember(target_desc,pcl_desc);
    %ss_same_moa.mat(ii,:) = contains(target_desc,pcl_desc);
    ss_same_moa.mat(ii,:) = ismember(target_desc,pcl_desc) | (contains(target_desc, '|') & contains(target_desc,pcl_desc));
    
    %check_matches = target_desc(~ismember(target_desc,pcl_desc) & contains(target_desc,pcl_desc));
    check_matches = target_desc(~ismember(target_desc,pcl_desc) & (contains(target_desc, '|') & contains(target_desc,pcl_desc)));
    
    if print_multi_target & length(check_matches) > 0
        pcl_desc
        check_matches
    end
end

sum(ss_same_moa.mat,'all')

ss_tbl.pcl_and_moa_agree = ss_same_moa.mat(:);




% Save the output table
%if ~isempty(project_id)
%    wtable(out_tbl,fullfile(wkdir,[project_id, '_correlation_to_selected_pcls.txt']));
%else
%    wtable(out_tbl,fullfile(wkdir,'correlation_to_selected_pcls.txt'));
%end

