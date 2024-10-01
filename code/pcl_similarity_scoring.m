function pcl_similarity_scoring(clusters_path,c_path,c_rank_path,col_meta_path,col_meta_kabx_path,outdir,min_clust_size,print_multi_target,stringify_cids)

% PCL_SIMILARITY_SCORING(CLUSTERS_PATH,DS_CORR,DS_CORR_RANK,COL_META,COL_META_KABX,OUTDIR,MIN_CLUST_SIZE,PRINT_MULTI_TARGET,STRINGIFY_CIDS)
% Set initial variables
% 
% datadir = '/idi/cgtb/morzech/idmp/screen4wk_screen5wk_kabx2_tbda1/analysis/pcls';
% gmtdir = fullfile(datadir, prefix, 'iteration_final');
% outdir = fullfile(datadir, prefix, 'iteration_final',['pcls_for_',project_id,'_final4']);
mk_cd_dir(outdir, false)

% Load clusters from spectral clustering
if ischar(clusters_path)
	clusters = parse_gmt(clusters_path);
end

% Set minimum cluster size if not provided in the input
if isempty(min_clust_size)
	min_clust_size = 2; 
end

% Set stringify_cids if not provided in the input
if isempty(stringify_cids)
	stringify_cids = false; 
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

if any(strcmp('pcl_desc', col_meta.Properties.VariableNames))
    col_meta = removevars(col_meta, 'pcl_desc');
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
%if isstr(c_rank_path)
%	c_rank = parse_gctx(c_rank_path);
%end
%c_rank = annotate_ds(c_rank,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','column','keyfield','cid');
%c_rank = annotate_ds(c_rank,table2struct(col_meta(:,{'cid','x_median_total_count'})),'dim','row','keyfield','cid');

% Make sure that c and c_rank are sorted in the same way
%c_rank = ds_slice(c_rank,'cid',c.cid,'rid',c.rid);

assert(all(diag(c.mat) == 1))

%assert(all(diag(c_rank.mat) == 1))

% Set the diagonal of the matrix to all NaN values
% - this will exclude self-similarity when calculating the median correlation of a treatment
% to a PCL that it is a member of
c.mat(eye(size(c.mat))==1) = NaN;

%c_rank.mat(eye(size(c_rank.mat))==1) = NaN;

assert(all(isnan(diag(c.mat))))

%assert(all(isnan(diag(c_rank.mat))))

% Add pert_id
% Already includes pert_id
%c.rdesc(:, c.rdict('pert_id')) = c.rdesc(:, c.rdict('broad_id'));
%idx = find_match(c.rdesc(:, c.rdict('broad_id')),'BRD-',true);
%c.rdesc(idx, c.rdict('pert_id')) = cellfun(@(s) s(1:13), c.rdesc(idx, c.rdict('broad_id')), 'uni', 0);

%% Check if all reference compounds are in the correlation matrix
setdiff(cat(1,clusters.entry),c.cid)

% Slice c to contain all profiles in rows and profiles from clusters in columns
ridx = ismember(c.rid, rids);
sum(ridx)
cidx = ismember(c.cid, cat(1,clusters.entry));
sum(cidx)

c = ds_slice(c,'ridx',ridx,'cidx',cidx);
%c_rank = ds_slice(c_rank,'ridx',ridx,'cidx',cidx);

% Save gctx files
disp('Saving gctx files')
mkgctx(fullfile(outdir,'ds_corr.gctx'), c)
%mkgctx(fullfile(outdir,'ds_corr_rank.gctx'), c_rank)

% ## Run MoA discovery using MoA-based clusters with overlap
%assert(isequal(c.cid, c_rank.cid),'Columns are not sorted identically')
%assert(isequal(c.rid, c_rank.rid),'Rows are not sorted identically')

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
    
    tmp_tbl.cid = c.rid;
    tmp_tbl.proj_broad_id = c.rdesc(:, c.rdict('proj_broad_id'));
    tmp_tbl.broad_id = c.rdesc(:, c.rdict('broad_id'));
    tmp_tbl.pert_id = c.rdesc(:, c.rdict('pert_id'));
    tmp_tbl.pert_idose = c.rdesc(:, c.rdict('pert_idose'));
    tmp_tbl.pert_dose = c.rdesc(:, c.rdict('pert_dose'));
    %tmp_tbl.x_median_total_count = c.rdesc(:, c.rdict('x_median_total_count'));

    num_cids = numel(c.rid);

    %clusters(jj).head

    tmp_tbl.cluster_id = repmat({clusters(jj).head}, num_cids, 1);

    tmp_tbl.cluster_desc = repmat({clusters(jj).desc}, num_cids, 1);

    tmp_tbl.cluster_size = repmat(clusters(jj).len, num_cids, 1);

    cidx = ismember(c.cid,clusters(jj).entry);

    tmp_tbl.cluster_size_actual = repmat(sum(cidx), num_cids, 1);

    if stringify_cids
        tmp_tbl.cluster_cids = repmat({stringify(unique(c.cid(cidx)))}, num_cids, 1);
    end
    tmp_tbl.cluster_proj_broad_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('proj_broad_id'))))}, num_cids, 1);
    tmp_tbl.cluster_broad_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('broad_id'))))}, num_cids, 1);
    tmp_tbl.cluster_pert_id = repmat({stringify(unique(c.cdesc(cidx,c.cdict('pert_id'))))}, num_cids, 1);
    
    tmp_tbl.cluster_num_pert_ids_unique = repmat(numel(unique(c.cdesc(cidx, c.cdict('pert_id')))), num_cids, 1);
    
    tmp_tbl.median_corr = median(c.mat(:, cidx), 2, 'omitnan');
    %tmp_tbl.median_corr_rank = median(c_rank.mat(:,cidx), 2, 'omitnan');

    tmp_tbl = struct2table(tmp_tbl);

    %tmp_tbl = sortrows(tmp_tbl,{'broad_id','median_corr_rank'},{'ascend','ascend'});
    
    out_tbl = [out_tbl;tmp_tbl];
    
end

size(out_tbl)

disp('All clusters successfully looped through')

headt(out_tbl)

% Make gctx
out_tbl.cid_idx = grp2idx(out_tbl.cid);
out_tbl.cluster_id_idx = grp2idx(out_tbl.cluster_id);

% Find indices
[~,ridx] = unique(out_tbl.cluster_id_idx);
rid = out_tbl.cluster_id(ridx);
[~,cidx] = unique(out_tbl.cid_idx);
cid = out_tbl.cid(cidx);

% Create metadata

%col_meta = out_tbl(cidx,{'cid','proj_broad_id','broad_id','pert_id','pert_idose','pert_dose'});
col_meta = col_meta;

if stringify_cids
    row_meta = out_tbl(ridx,{'cluster_id','cluster_desc','cluster_size','cluster_cids','cluster_proj_broad_id','cluster_broad_id','cluster_pert_id','cluster_num_pert_ids_unique'});
else
    row_meta = out_tbl(ridx,{'cluster_id','cluster_desc','cluster_size','cluster_proj_broad_id','cluster_broad_id','cluster_pert_id','cluster_num_pert_ids_unique'});
end

% Create an empty matrix
mat = nan(numel(rid),numel(cid));
ind = sub2ind(size(mat),out_tbl.cluster_id_idx,out_tbl.cid_idx);

ds = mkgctstruct(mat,'rid',rid','cid',cid);
ds = annotate_ds(ds,table2struct(col_meta),'dim','column','keyfield','cid');
ds = annotate_ds(ds,table2struct(row_meta),'dim','row','keyfield','cluster_id');

% Create gctx files with various data
%fields = {'median_corr','median_corr_rank'};
fields = {'median_corr'};

for ii = 1:numel(fields)
    disp(fields{ii})
    ds_tmp = ds;
    ds_tmp.mat(ind) = out_tbl.(fields{ii});

    % Verify
    rid1 = ds_tmp.rid(2);
    cid1 = ds_tmp.cid(1);
    corr1 = ds_tmp.mat(2,1);
    idx = ismember(out_tbl.cid,cid1)&ismember(out_tbl.cluster_id,rid1);
    assert(sum(idx)>0,'Rows or columns of the gct structure do not match');
    assert(out_tbl.(fields{ii})(idx)==corr1,'Value in the matrix do not match the one in the table');
    mkgctx(fullfile(outdir,['cluster_',fields{ii},'.gctx']), ds_tmp)
end

% Section to check clusters for their most similar KABX dsCGI profile and define them as PCLs if the profiles are in-MOA (share same MOA as PCL)
% and exclude them as uninterpretable clusters for MOA prediction otherwise
disp('Beginning check on clusters for their most similar KABX dsCGI profile and define them as PCLs if the profiles are in-MOA (share same MOA as PCL) and exclude them as uninterpretable clusters for MOA prediction otherwise');

% Parse gctx with the cluster similarity score
g = glob(fullfile(outdir,'cluster_median_corr_n*.gctx'))
ss_all = parse_gctx(g{1}); % cluster similarity score of all treatments

col_meta_ss_all = cell2table([ss_all.cid,ss_all.cdesc],'VariableNames',['cid';ss_all.chd]);
col_meta_ss_all.target_description = any2str(col_meta_ss_all.target_description);
disp('Head of col_meta_ss_all tbl:');
disp(headt(col_meta_ss_all));
disp('Unique KABX MOA annotations including multi-target:');
unique(col_meta_ss_all.target_description)

row_meta_ss_all = cell2table([ss_all.rid,ss_all.rdesc],'VariableNames',['rid';ss_all.rhd]);
row_meta_ss_all.cluster_desc = any2str(row_meta_ss_all.cluster_desc);
disp('Head of row_meta_ss_all tbl:');
disp(headt(row_meta_ss_all));
disp('Unique cluster MOA annotations:');
unique(row_meta_ss_all.cluster_desc)

[a,b] = ind2sub(size(ss_all.mat),1:numel(ss_all.mat));

ss_all_tbl = [col_meta_ss_all(b,:),row_meta_ss_all(a,:)];
disp('Size of ss_all tbl:');
disp(size(ss_all_tbl));
ss_all_tbl.cluster_similarity_score = ss_all.mat(:);

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

disp('Sum of treatment-PCL pairs with matching MOA (binary label of 1):');
disp(sum(ss_all_same_moa.mat,'all'));

ss_all_tbl.cluster_and_moa_agree = ss_all_same_moa.mat(:);

% Load KABX column metadata
if isstr(col_meta_kabx_path)
	col_meta_kabx = rtable(col_meta_kabx_path);
end

% Select profiles with annotated MOA from KABX
% cidx = ~cellfun(@isempty,s.cdesc(:,s.cdict('target_description')));
%cidx = ismember(s.cid, col_meta_kabx.cid);
cidx = ismember(ss_all_tbl.cid, col_meta_kabx.cid);
disp('Number of rows of ss_all_tbl from KABX treatments:');
disp(sum(cidx));

%ss_kabx = ds_slice(ss_all,'cidx',cidx)
ss_kabx_tbl = ss_all_tbl(cidx,:);

disp('Size of ss_all tbl:');
disp(size(ss_all_tbl));

disp(sprintf('Number of pert_ids in ss_all_tbl: %d', numel(unique(ss_all_tbl.pert_id))));
disp(sprintf('Number of broad_ids in ss_all_tbl: %d', numel(unique(ss_all_tbl.broad_id))));
disp(sprintf('Number of proj_broad_ids in ss_all_tbl: %d', numel(unique(ss_all_tbl.proj_broad_id))));
disp(sprintf('Number of treatments (cids) in ss_all_tbl: %d', numel(unique(ss_all_tbl.cid))));
disp(sprintf('Number of treatment-PCL pairs with matching MOA (binary label of 1) in ss_all_tbl: %d', sum(ss_all_tbl.cluster_and_moa_agree)));

disp('Size of ss_kabx tbl:');
disp(size(ss_kabx_tbl));

disp(sprintf('Number of pert_ids in ss_kabx_tbl: %d', numel(unique(ss_kabx_tbl.pert_id))));
disp(sprintf('Number of broad_ids in ss_kabx_tbl: %d', numel(unique(ss_kabx_tbl.broad_id))));
disp(sprintf('Number of proj_broad_ids in ss_kabx_tbl: %d', numel(unique(ss_kabx_tbl.proj_broad_id))));
disp(sprintf('Number of treatments (cids) in ss_kabx_tbl: %d', numel(unique(ss_kabx_tbl.cid))));
disp(sprintf('Number of treatment-PCL pairs with matching MOA (binary label of 1) in ss_kabx_tbl: %d', sum(ss_kabx_tbl.cluster_and_moa_agree)));


cluster_list = unique(ss_kabx_tbl.rid);
disp(sprintf('Number of clusters: %d', length(cluster_list)));

moa_cluster_list = unique(ss_kabx_tbl.cluster_desc);
disp(sprintf('Number of MOAs represented in clusters: %d', length(moa_cluster_list)));

% Filter out uninterpretable clusters (non MOA-separable) to find PCL clusters whose most similar treatments are in-MOA rather than out-of-MOA
disp('Using ss_kabx_tbl to find most similar KABX treatment to each cluster');

% Find the row with the maximum cluster_similarity_score for each rid

maxSimilarityScore = groupsummary(ss_kabx_tbl, 'rid', @max, 'cluster_similarity_score');

size(maxSimilarityScore)

head(maxSimilarityScore)

disp('Size of ss_kabx tbl:');
disp(size(ss_kabx_tbl));

ss_kabx_tbl = innerjoin(ss_kabx_tbl, maxSimilarityScore, 'Keys', 'rid');

disp('Size of ss_kabx tbl:');
disp(size(ss_kabx_tbl));

disp('Head of ss_kabx tbl:');
disp(headt(ss_kabx_tbl));

% Get the rows with the maximum cluster_similarity_score
maxRows = ss_kabx_tbl(ss_kabx_tbl.cluster_similarity_score >= ss_kabx_tbl.fun1_cluster_similarity_score,:);

size(maxRows)

headt(maxRows)

% Find the rows where cluster_and_moa_agree = 1
pclListRows = maxRows(maxRows.cluster_and_moa_agree == 1,:);

% Get the rids of these rows
pcl_list = unique(pclListRows.rid);

disp(sprintf('Length of PCL cluster list: %d', length(pcl_list)));

moa_pcl_list = unique(pclListRows.cluster_desc);
disp(sprintf('Number of MOAs represented in PCL cluster list: %d', length(moa_pcl_list)));

% Remove clusters whose most similar KABX dsCGI profile is of a different MOA reflecting signal not specific to one particular MOA
% while clusters whose most similar KABX dsCGI profile is in-MOA are defined as PCL clusters and used for making MOA predictions
disp('Removing clusters whose most similar KABX dsCGI profile is of a different MOA reflecting signal not specific to one particular MOA');
disp('Clusters whose most similar KABX dsCGI profile is in-MOA are defined as PCL clusters and will be used for making MOA predictions');

num_clusters_ini = length(cluster_list);
num_pcl_clusters = length(pcl_list);
num_uninterpretable_clusters = num_clusters_ini - num_pcl_clusters;
disp(sprintf('The initial number of clusters is %d', num_clusters_ini));
disp(sprintf('Removing %d uninterpretable clusters with most similar KABX dsCGI profile out of MOA', num_uninterpretable_clusters));
disp(sprintf('The remaining number of PCL clusters is %d', num_pcl_clusters));

pcls = parse_gmt(clusters_path);
pcls_head = struct2table(pcls).head;
length(pcls_head)
pcls(~ismember(pcls_head, pcl_list)) = []; % remove uninterpretable clusters 
pcls_tbl = struct2table(pcls);
disp('Head and size of PCLs tbl');
disp(headt(pcls_tbl));
disp(size(pcls_tbl));

mkgmt(fullfile(outdir, 'pcls.gmt'), pcls)

disp('Size of ss_all tbl:');
disp(size(ss_all_tbl));

ss_all_tbl = innerjoin(ss_all_tbl, maxSimilarityScore, 'Keys', 'rid');

ss_all_tbl = ss_all_tbl(ismember(ss_all_tbl.rid, pcl_list),:);

disp('Size of ss_all tbl:');
disp(size(ss_all_tbl));

if any(strcmp('pcl_desc', ss_all_tbl.Properties.VariableNames))
    ss_all_tbl = removevars(ss_all_tbl, 'pcl_desc');
end

if any(strcmp('GroupCount', ss_all_tbl.Properties.VariableNames))
    ss_all_tbl = removevars(ss_all_tbl, 'GroupCount');
end

if any(strcmp('fun1_cluster_similarity_score', ss_all_tbl.Properties.VariableNames))
    ss_all_tbl.cluster_max_kabx_similarity_score =  ss_all_tbl.fun1_cluster_similarity_score;

    ss_all_tbl = removevars(ss_all_tbl, 'fun1_cluster_similarity_score');
end

% Get the current column names of the table
col_names = ss_all_tbl.Properties.VariableNames;

% Replace "cluster" with "pcl" in each column name
new_col_names = strrep(col_names, 'cluster', 'pcl');

% Assign the updated column names back to the table
ss_all_tbl.Properties.VariableNames = new_col_names;

ss_all_tbl.pcl_id = ss_all_tbl.rid;

disp('Size and head of ss_all tbl:');
disp(size(ss_all_tbl));
disp(headt(ss_all_tbl));

% Make gctx
ss_all_tbl.cid_idx = grp2idx(ss_all_tbl.cid);
ss_all_tbl.pcl_id_idx = grp2idx(ss_all_tbl.pcl_id);

% Find indices
[~,ridx] = unique(ss_all_tbl.pcl_id_idx);
rid = ss_all_tbl.pcl_id(ridx);
[~,cidx] = unique(ss_all_tbl.cid_idx);
cid = ss_all_tbl.cid(cidx);

% Create metadata

col_names = ss_all_tbl.Properties.VariableNames;

% Find all column names that do not contain "pcl"
non_pcl_col_names = col_names(~contains(col_names, 'pcl') & ~contains(col_names, {'rid','cid_idx'}));
pcl_col_names = col_names(contains(col_names, 'pcl') & ~contains(col_names, {'pcl_similarity_score','pcl_and_moa_agree','pcl_id_idx'}));

col_meta = ss_all_tbl(cidx, non_pcl_col_names);

disp('Size and head of final col_meta for gct files:');
disp(size(col_meta));
disp(headt(col_meta));

row_meta = ss_all_tbl(ridx, pcl_col_names);

disp('Size and head of final row_meta for gct files:');
disp(size(row_meta));
disp(headt(row_meta));

% Create an empty matrix
mat = nan(numel(rid),numel(cid));
ind = sub2ind(size(mat),ss_all_tbl.pcl_id_idx,ss_all_tbl.cid_idx);

ds = mkgctstruct(mat,'rid',rid','cid',cid);
ds = annotate_ds(ds,table2struct(col_meta),'dim','column','keyfield','cid');
ds = annotate_ds(ds,table2struct(row_meta),'dim','row','keyfield','pcl_id');

% Create gctx files with various data
fields = {'pcl_similarity_score','pcl_and_moa_agree'};

for ii = 1:numel(fields)
    disp(fields{ii})
    ds_tmp = ds;
    ds_tmp.mat(ind) = ss_all_tbl.(fields{ii});

    % Verify
    rid1 = ds_tmp.rid(2);
    cid1 = ds_tmp.cid(1);
    corr1 = ds_tmp.mat(2,1);
    idx = ismember(ss_all_tbl.cid,cid1)&ismember(ss_all_tbl.pcl_id,rid1);
    assert(sum(idx)>0,'Rows or columns of the gct structure do not match');
    assert(ss_all_tbl.(fields{ii})(idx)==corr1,'Value in the matrix do not match the one in the table');
    mkgctx(fullfile(outdir,[fields{ii},'.gctx']), ds_tmp)
end

% Save the output table
%if ~isempty(project_id)
%    wtable(ss_all_tbl,fullfile(outdir,[project_id, '_correlation_to_selected_pcls.txt']));
%else
%    wtable(ss_all_tbl,fullfile(outdir,'correlation_to_selected_pcls.txt'));
%end

