function pcl_similarity_scoring(clusters_path,c_path,c_rank_path,col_meta_path,fgr_path,wkdir,min_clust_size, col_meta_kabx_path, print_multi_target)

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
% - this will exclude self-similarity when calculating the median correlation of a treatment
% to a PCL that it is a member of
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
    mkgctx(fullfile(wkdir,['cluster_',fields{ii},'.gctx']), ds_tmp)
end

% Section to check clusters for their most similar KABX dsCGI profile and define them as PCLs if the profiles are in-MOA (share same MOA as PCL)
% and exclude them as uninterpretable clusters for MOA prediction otherwise

% Parse gctx with the cluster similarity score
g = glob(fullfile(wkdir,'cluster_median_corr_n*.gctx'))
ss_all = parse_gctx(g{1}); % cluster similarity score of all treatments

col_meta_ss_all = cell2table([ss_all.cid,ss_all.cdesc],'VariableNames',['cid';ss_all.chd]);
col_meta_ss_all.target_description = any2str(col_meta_ss_all.target_description);
headt(col_meta_ss_all)
unique(col_meta_ss_all.target_description)

row_meta_ss_all = cell2table([ss_all.rid,ss_all.rdesc],'VariableNames',['rid';ss_all.rhd]);
row_meta_ss_all.cluster_desc = any2str(row_meta_ss_all.cluster_desc);
headt(row_meta_ss_all)
unique(row_meta_ss_all.cluster_desc)

[a,b] = ind2sub(size(ss_all.mat),1:numel(ss_all.mat));

ss_all_tbl = [col_meta_ss_all(b,:),row_meta_ss_all(a,:)];
size(ss_all_tbl)
ss_all_tbl.similarity_score = ss_all.mat(:);

ss_all_same_moa = ss_all;
ss_all_same_moa.mat = zeros(size(ss_all.mat));
target_desc = ss_all_same_moa.cdesc(:, ss_all_same_moa.cdict('target_description'));

% Set print_multi_target if not provided in the input
if isempty(print_multi_target)
	print_multi_target = false; 
end

for ii = 1:numel(ss_all_same_moa.rid)
    cluster_desc = ss_all_same_moa.rdesc(ii,ss_all_same_moa.rdict('cluster_desc'));
    %ss_all_same_moa.mat(ii,:) = ismember(target_desc,cluster_desc);
    %ss_all_same_moa.mat(ii,:) = contains(target_desc,cluster_desc);
    ss_all_same_moa.mat(ii,:) = ismember(target_desc,cluster_desc) | (contains(target_desc, '|') & contains(target_desc,cluster_desc));
    
    %check_matches = target_desc(~ismember(target_desc,cluster_desc) & contains(target_desc,cluster_desc));
    check_matches = target_desc(~ismember(target_desc,cluster_desc) & (contains(target_desc, '|') & contains(target_desc,cluster_desc)));
    
    if print_multi_target & length(check_matches) > 0
        cluster_desc
        check_matches
    end
end

sum(ss_all_same_moa.mat,'all')

ss_all_tbl.cluster_and_moa_agree = ss_all_same_moa.mat(:);


% Load KABX column metadata
if isstr(col_meta_kabx_path)
	col_meta_kabx = rtable(col_meta_kabx_path);
end

% Select profiles with annotated MOA from KABX
% cidx = ~cellfun(@isempty,s.cdesc(:,s.cdict('target_description')));
%cidx = ismember(s.cid, col_meta_kabx.cid);
cidx = ismember(ss_all_tbl.cid, col_meta_kabx.cid);
sum(cidx)

%ss_kabx = ds_slice(ss_all,'cidx',cidx)
ss_kabx_tbl = ss_all_tbl(cidx,:);

size(ss_all_tbl)

numel(unique(ss_all_tbl.pert_id))
numel(unique(ss_all_tbl.broad_id))
numel(unique(ss_all_tbl.proj_broad_id))
numel(unique(ss_all_tbl.cid))
sum(ss_all_tbl.cluster_and_moa_agree)

size(ss_kabx_tbl)

numel(unique(ss_kabx_tbl.pert_id))
numel(unique(ss_kabx_tbl.broad_id))
numel(unique(ss_kabx_tbl.proj_broad_id))
numel(unique(ss_kabx_tbl.cid))
sum(ss_kabx_tbl.cluster_and_moa_agree)


cluster_list = unique(ss_kabx_tbl.rid);
length(cluster_list)

moa_cluster_list = unique(ss_kabx_tbl.cluster_desc);
length(moa_cluster_list)

% filter out MOA in-separable PCLs prior to iterating

% Find the row with the maximum similarity_score for each rid

maxSimilarityScore = groupsummary(ss_kabx_tbl, 'rid', @max, 'similarity_score');

size(maxSimilarityScore)

head(maxSimilarityScore)

size(ss_kabx_tbl)

ss_kabx_tbl = innerjoin(ss_kabx_tbl, maxSimilarityScore, 'Keys', 'rid');

size(ss_kabx_tbl)

headt(ss_kabx_tbl)

% Get the rows with the maximum similarity_score
maxRows = ss_kabx_tbl(ss_kabx_tbl.similarity_score >= ss_kabx_tbl.fun1_similarity_score,:);

size(maxRows)

headt(maxRows)

% Find the rows where cluster_and_moa_agree = 1
pclListRows = maxRows(maxRows.cluster_and_moa_agree == 1,:);

% Get the rids of these rows
pcl_list = unique(pclListRows.rid);

length(pcl_list)

moa_pcl_list = unique(pclListRows.cluster_desc);
length(moa_pcl_list)

% Remove clusters whose most similar KABX dsCGI profile is of a different MOA reflecting signal not specific to one particular MOA
% while clusters whose most similar KABX dsCGI profile is in-MOA are defined as PCL clusters and used for making MOA predictions
num_clusters_ini = length(cluster_list);
num_pcl_clusters = length(pcl_list);
num_uninterpretable_clusters = num_clusters_ini - num_pcl_clusters;
disp(sprintf('The initial number of clusters is %d', num_clusters_ini))
disp(sprintf('Removing %d uninterpretable clusters with most similar KABX dsCGI profile out of MOA', num_uninterpretable_clusters))
disp(sprintf('The remaining number of PCL clusters is %d', num_pcl_clusters))

pcls = parse_gmt(clusters_path);
pcls(~ismember(pcls.head, pcl_list)) = []; % remove uninterpretable clusters 
pcls_tbl = struct2table(pcls);
headt(pcls_tbl)
size(pcls_tbl)

mkgmt(fullfile(wkdir, 'pcls.gmt'), pcls)




% Save the output table
%if ~isempty(project_id)
%    wtable(out_tbl,fullfile(wkdir,[project_id, '_correlation_to_selected_pcls.txt']));
%else
%    wtable(out_tbl,fullfile(wkdir,'correlation_to_selected_pcls.txt'));
%end

