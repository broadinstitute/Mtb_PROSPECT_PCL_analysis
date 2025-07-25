function pcl_confidence_scoring(pcls_path,similarity_score_path,pcl_and_moa_agree_labels_path,col_meta_reference_set_for_pcls_path,high_confidence_pcl_confidence_score_thres,outdir,unknown_target_description_values,make_fig,out_tbl_savename,opt_tbl_savename,out_tbl_test_cmpd_savename)

%% Check inputs %%

mk_cd_dir(outdir, false);

%exist(pcls_path) > 0

%exist(similarity_score_path) > 0

%exist(pcl_and_moa_agree_labels_path) > 0

%exist(col_meta_reference_set_for_pcls_path) > 0

assert(exist(pcls_path) > 0, 'PCLs path does not exist')

assert(exist(similarity_score_path) > 0, 'PCL similarity scores path does not exist')

assert(exist(pcl_and_moa_agree_labels_path) > 0, 'pcl_and_moa_agree path does not exist')

assert(exist(col_meta_reference_set_for_pcls_path) > 0, 'col_meta_reference_set_for_pcls path does not exist')

%ls(outdir)

% Set high_confidence_pcl_confidence_score_thres if not provided in the input
if isempty(high_confidence_pcl_confidence_score_thres)
	high_confidence_pcl_confidence_score_thres = 1; % defined as such for uniformity and simplicity across all PCL clusters
end

% Set unknown_target_description_values if not provided in the input
if isempty(unknown_target_description_values)
	unknown_target_description_values = {'','NA','NAN','NaN','whole cell only','unknown'}; % in case of different processing filetypes and NA values stored in import/export
end

% Set make_fig if not provided in the input
if isempty(make_fig)
	make_fig = false;
end

% Set out_tbl_savename if not provided in the input
if isempty(out_tbl_savename)
	out_tbl_savename = 'by_pcl_similarity_to_confidence_score_thresholds_from_training_on_reference_set.txt';
end

% Set opt_tbl_savename if not provided in the input
if isempty(opt_tbl_savename)
	opt_tbl_savename = 'by_pcl_high_confidence_similarity_score_thresholds_from_training_on_reference_set.txt';
end

% Set out_tbl_test_cmpd_savename if not provided in the input
if isempty(out_tbl_test_cmpd_savename)
	out_tbl_test_cmpd_savename = 'by_pcl_similarity_to_confidence_score_test_cmpd_results.txt';
end


% ## Load annotation of reference set CGI profiles used to construct PCL clusters and extract the unique reference set compound list (pert_ids)

%% will be used to split up training (reference set) set from test (unknown compound) set and use of col_meta_reference_set_for_pcls.txt file
%% is particularly important for iterations of LOOCV in order to subset out the leftout reference set compound as the test set since
%% it was not used in PCL cluster construction and will not appear in the extracted reference set compound list 

col_meta_reference_set_for_pcls_path
col_meta_reference_set_for_pcls = rtable(col_meta_reference_set_for_pcls_path);
disp('Size and head of col_meta_reference_set_for_pcls table');
size(col_meta_reference_set_for_pcls)
headt(col_meta_reference_set_for_pcls)

unique_pert_ids_in_reference_set_for_pcls = unique(col_meta_reference_set_for_pcls.pert_id);

disp(sprintf('Number of compounds (pert_ids) in reference set used for PCL cluster construction: %d', length(unique_pert_ids_in_reference_set_for_pcls)))

% ## Parse gctx with the PCL cluster similarity score (median correlation of condition to each PCL cluster's conditions)

%g = glob(fullfile(outdir,'pcl_similarity_score_n*.gctx'))
%s = parse_gctx(g{1});
s = parse_gctx(similarity_score_path);

s.chd

% ## Parse gctx with binary labels if reference set condition belongs to MOA of a PCL cluster, i.e, pcl_and_moa_agree

%g = glob(fullfile(outdir,'pcl_and_moa_agree_n*.gctx'))
%l = parse_gctx(g{1});
l = parse_gctx(pcl_and_moa_agree_labels_path);

l.chd

% ## Select all profiles regardless of rcorr

disp(sprintf('Total number of compounds (pert_ids): %d', numel(unique(s.cdesc(:,s.cdict('pert_id'))))))
disp(sprintf('Total number of broad_ids: %d', numel(unique(s.cdesc(:,s.cdict('broad_id'))))))
disp(sprintf('Total number of proj_broad_ids: %d', numel(unique(s.cdesc(:,s.cdict('proj_broad_id'))))))
disp(sprintf('Total number of CGI profiles / conditions (cids): %d', numel(unique(s.cid))))

%% Commented out portion would consider only CGI profiles with rcorr (replicate correlation/reproducibility) above a given threshold %%
%cidx = cell2mat(s.cdesc(:,s.cdict('rcorr')))>=0.3;
%sum(cidx)
%disp(sprintf('Total number of compounds (pert_ids) with high rcorr: %d', numel(unique(s.cdesc(cidx,s.cdict('pert_id'))))))
%disp(sprintf('Total number of broad_ids with high rcorr: %d', numel(unique(s.cdesc(cidx,s.cdict('broad_id'))))))
%disp(sprintf('Total number of proj_broad_ids with high rcorr: %d', numel(unique(s.cdesc(cidx,s.cdict('proj_broad_id'))))))
%disp(sprintf('Total number of CGI profiles / conditions (cids) with high rcorr: %d', numel(unique(s.cdesc(cidx,s.cdict('pert_id'))))))
%
%s = ds_slice(s,'cidx',cidx);

unknown_target_description_values

cidx = ~cellfun(@(x) any(strcmp(x, unknown_target_description_values)), s.cdesc(:, s.cdict('target_description')));
disp('Number of MOA-annotated CGI profiles')
sum(cidx)
disp('Number of CGI profiles lacking MOA annotation')
sum(~cidx)

% CGI profiles from MOA-annotated reference set compounds
sortrows(cell2table(tabulate(s.cdesc(cidx,s.cdict('target_description'))),...
    'VariableNames',{'target_description','count','percent'}),1)

% CGI profiles from unannotated test compounds NA, N/A, NAN, Not Applicable

disp('Unique target_description values for unknown conditions:');
unique(s.cdesc(~cidx,s.cdict('target_description')))
disp('Setting unknown target_description values to NA');
s.cdesc(~cidx,s.cdict('target_description')) = {'NA'};
disp('Now unique target_description values for unknown conditions:');
unique(s.cdesc(~cidx,s.cdict('target_description')))

sortrows(cell2table(tabulate(s.cdesc(~cidx,s.cdict('target_description'))),...
    'VariableNames',{'target_description','count','percent'}),1)

% ## Construct table with PCL cluster similarity score and binary labels of if CGI profile matches each PCL cluster's annotated MOA for fitting PCL cluster confidence score mapping and finding high-confidence similarity score thresholds

ss = s

% Get column and row metadata fields for table
col_meta_ss = cell2table([ss.cid,ss.cdesc],'VariableNames',['cid';ss.chd]);
col_meta_ss.target_description = any2str(col_meta_ss.target_description);
disp('Head of col_meta_ss tbl:');
disp(headt(col_meta_ss));
disp('Unique reference set MOA annotations including multi-target:');
unique(col_meta_ss.target_description)

row_meta_ss = cell2table([ss.rid,ss.rdesc],'VariableNames',['rid';ss.rhd]);
row_meta_ss.pcl_desc = any2str(row_meta_ss.pcl_desc);
disp('Head of row_meta_ss tbl:');
disp(headt(row_meta_ss));
disp('Unique PCL cluster MOA annotations:');
unique(row_meta_ss.pcl_desc)

[a,b] = ind2sub(size(ss.mat),1:numel(ss.mat));

ss_tbl = [col_meta_ss(b,:),row_meta_ss(a,:)];
disp('Size of ss_tbl with column and row metadata fields from similarity score gct');
disp(size(ss_tbl));
disp('Attaching PCL cluster similarity score');
ss_tbl.pcl_similarity_score = ss.mat(:);
disp(size(ss_tbl));

% %% Commented out code section to set target_description values for unknown compounds to empty cell as opposed to NA, NAN, or other unknown_target_description_values %%
% unknown_target_description_values
% idx_unknown = ismember(ss_tbl.target_description,unknown_target_description_values);
% disp('Number of CGI profile rows lacking MOA annotation:');
% sum(idx_unknown)
% disp('Unique target_description values for unknown conditions:');
% unique(ss_tbl.target_description(idx_unknown))
% disp('Setting unknown target_description values to empty cell');
% ss_tbl.target_description(idx_unknown) = {''};
% disp('Now unique target_description values for unknown conditions:');
% unique(ss_tbl.target_description(idx_unknown))
% disp(size(ss_tbl));

% ### Add binary CGI profile-PCL MOA agreement labels from pcl_and_moa_agree gct

l = ds_slice(l, 'cid', ss.cid, 'rid', ss.rid);

assert(isequal(ss.rid,l.rid), 'ss and l rows not in matching order')
assert(isequal(ss.cid,l.cid), 'ss and l rows not in matching order')

disp('Attaching CGI profile-PCL MOA agreement binary labels (1 for in-MOA / 0 for out-of-MOA)');
ss_tbl.pcl_and_moa_agree = l.mat(:);
disp(size(ss_tbl));

% ### Rank PCL clusters (highest to lowest PCL cluster similarity score) for every CGI profile / condition

ss_rank = ss;
ss_rank.mat = rankorder(ss.mat,'descend');

assert(isequal(ss.rid,ss_rank.rid), 'ss and ss_rank rows not in matching order')
assert(isequal(ss.cid,ss_rank.cid), 'ss and ss_rank rows not in matching order')

disp('Attaching rank of PCL clusters for each CGI profile from high-to-low PCL cluster similarity score');
ss_tbl.rank_pcl = ss_rank.mat(:);
disp(size(ss_tbl));

% ### Identify any MOAs that do not have matching PCL cluster

all_targets = unique(ss_tbl.target_description);
length(all_targets)

targets_with_match = unique(ss_tbl.target_description(ss_tbl.pcl_and_moa_agree == 1));
length(targets_with_match)

targets_without_match = all_targets(~ismember(all_targets, targets_with_match))
length(targets_without_match)

% ## Split training and test set

ss_tbl_train_and_test = ss_tbl;

disp('Size of train and test set tbl:');
disp(size(ss_tbl_train_and_test));

test_cmpd_idx = ~ismember(ss_tbl_train_and_test.pert_id, unique_pert_ids_in_reference_set_for_pcls);

disp(sprintf('Number of rows that are test (unknown compound) conditions: %d', sum(test_cmpd_idx)));

ss_tbl_test_cmpd = ss_tbl_train_and_test(test_cmpd_idx, :);

disp('Size and head of test set tbl:');
disp(size(ss_tbl_test_cmpd));
disp(headt(ss_tbl_test_cmpd));

disp(sprintf('Number of rows that are training (reference set) conditions: %d', sum(~test_cmpd_idx)));

ss_tbl = ss_tbl_train_and_test(~test_cmpd_idx, :);

disp('Size and head of training (reference) set tbl:');
disp(size(ss_tbl));
disp(headt(ss_tbl));

assert(all(ismember(unique(ss_tbl.pert_id), unique_pert_ids_in_reference_set_for_pcls)), 'One or more compounds in ss_tbl are not from reference set annotation')

assert(all(ismember(unique_pert_ids_in_reference_set_for_pcls, unique(ss_tbl.pert_id))), 'One or more compounds from reference set annotation are not in ss_tbl')

assert(size(ss_tbl, 1) + size(ss_tbl_test_cmpd, 1) == size(ss_tbl_train_and_test, 1), 'Number of rows of ss_tbl and ss_tbl_test_cmpd do not equal total number of rows of ss_tbl_train_and_test')

% ## Run ROC analysis of PCL cluster similarity score to MOA agreement over all CGI profiles-PCL pairs

% Define colors: gray for 0 and purple for 1
if make_fig
    colors = [0.8 0.8 0.8; 0.5 0 0.5];
end

if make_fig
    figure
    gscatter(ss_tbl.pcl_similarity_score, ss_tbl.rcorr, ss_tbl.pcl_and_moa_agree, colors, '.', 8)
    xlabel('PCL similarity score')
    ylabel('Replicate correlation')
    title(sprintf('All reference set CGI profile-PCL pairs\nOvR binary classifier label: 1 for in-MOA / 0 for out-of-MOA'))
    %subtitle('OvR binary classifer label: 1 for in-MOA / 0 for out-of-MOA')
    %xlim([-1, 1])
    %ylim([-1, 1])

    saveas_png(gcf, outdir, ['reference_set_pcl_similarity_score_x_rcorr_all_profile_pcl_pairs.png'])
    close(gcf)
end

if make_fig
    [x,y,t,auc,optrocpt] = perfcurve(ss_tbl.pcl_and_moa_agree,ss_tbl.pcl_similarity_score,1,'NegClass',0,'UseNearest','off','XVals','all','Cost',[0 1; 1 0]);
    auc
    optrocpt

    figure
        plot(x, y)
        xlabel('False positive rate')
        ylabel('True positive rate')
        title(sprintf('ROC analysis of PCL cluster similarity score to\nPCL-MOA agreement for all reference set CGI profile-PCL pairs'))
        text(0.5, 0.5, ['AUC = ' num2str(auc, 3)], 'FontSize', 12, 'Color', 'black')
        text(0.5, 0.4, ['Opt ROC Pt = (' num2str(optrocpt(1), 3) ', ' num2str(optrocpt(2), 3) ')'], 'FontSize', 12, 'Color', 'black')
    saveas_png(gcf, outdir, ['reference_set_roc_fpr_x_tpr_all_profile_pcl_pairs.png'])
    close(gcf)
    
    figure
        plot(y,t)
        xlabel('True positive rate')
        ylabel('PCL similarity score threshold')
        title(sprintf('ROC analysis of PCL cluster similarity score to\nPCL-MOA agreement for all reference set CGI profile-PCL pairs'))
    saveas_png(gcf, outdir, ['reference_set_roc_tpr_x_pcl_similarity_score_threshold_all_profile_pcl_pairs.png'])
    close(gcf)
end

if make_fig
    figure
    gscatter(ss_tbl_test_cmpd.pcl_similarity_score,ss_tbl_test_cmpd.rcorr,ss_tbl_test_cmpd.pcl_and_moa_agree, colors, '.', 8)
    xlabel('PCL similarity score')
    ylabel('Replicate correlation')
    title(sprintf('All test compound CGI profile-PCL pairs\nAll 0 for unknown MOA'))
    %subtitle('0 for unknown MOA')
    %xlim([-1, 1])
    %ylim([-1, 1])

    saveas_png(gcf, outdir, ['test_compound_pcl_similarity_score_x_rcorr_all_profile_pcl_pairs.png'])
    close(gcf)
end

% ## Check for any missing PCL cluster similarity scores in the training (reference) set

idx = ~isnan(ss_tbl.pcl_similarity_score);
disp(sprintf('Total number of reference set CGI profile-PCL pairs with PCL cluster similarity score to be used in training: %d', sum(idx)));
ss_tbl_sele = ss_tbl(idx,:);

assert(sum(~idx) == 0, sprintf('There are reference set CGI profiles with missing PCL cluster similarity scores: %d rows', sum(~idx)))
assert(sum(idx) == height(ss_tbl), sprintf('There are reference set CGI profiles with missing PCL cluster similarity scores: %d rows', height(ss_tbl) - sum(idx)))

% ## Verify number of training (reference) set compounds and conditions

disp(sprintf('Number of pert_ids in training (reference) set ss_tbl: %d', numel(unique(ss_tbl_sele.pert_id))));
disp(sprintf('Number of broad_ids in training (reference) set ss_tbl: %d', numel(unique(ss_tbl_sele.broad_id))));
disp(sprintf('Number of proj_broad_ids in training (reference) set ss_tbl: %d', numel(unique(ss_tbl_sele.proj_broad_id))));
disp(sprintf('Number of conditions (cids) in training (reference) set ss_tbl: %d', numel(unique(ss_tbl_sele.cid))));

% ## Verify number of test (unknown compound) set compounds and conditions

disp(sprintf('Number of pert_ids in test set ss_tbl_test_cmpd: %d', numel(unique(ss_tbl_test_cmpd.pert_id))));
disp(sprintf('Number of broad_ids in test set ss_tbl_test_cmpd: %d', numel(unique(ss_tbl_test_cmpd.broad_id))));
disp(sprintf('Number of proj_broad_ids in test set ss_tbl_test_cmpd: %d', numel(unique(ss_tbl_test_cmpd.proj_broad_id))));
disp(sprintf('Number of conditions (cids) in test set ss_tbl_test_cmpd: %d', numel(unique(ss_tbl_test_cmpd.cid))));

% ## Verify number of PCL clusters in reference set and test set tables matches gmt file from previous steps

pcl_list = unique(ss_tbl_sele.rid);
disp(sprintf('Number of PCL clusters in reference set table ss_tbl_sele: %d', length(pcl_list)));

moa_list = unique(ss_tbl_sele.pcl_desc);
disp(sprintf('Number of MOAs represented by PCL clusters in reference set table ss_tbl_sele: %d', length(moa_list)));

pcls = parse_gmt(pcls_path);
pcls_tbl = struct2table(pcls);
disp('Head and size of pcls_tbl read in from PCL cluster gmt file:');
headt(pcls_tbl)
size(pcls_tbl)

assert(height(pcls_tbl) == length(pcl_list), 'Number of PCL clusters in gmt file and reference set table do not match')

disp(sprintf('Number of PCL clusters in test set table ss_tbl_test_cmpd: %d', length(unique(ss_tbl_test_cmpd.rid))));

disp(sprintf('Number of MOAs represented by PCL clusters in test set table ss_tbl_test_cmpd: %d', length(unique(ss_tbl_test_cmpd.pcl_desc))));

assert(all(ismember(ss_tbl_test_cmpd.rid, pcl_list)), 'Number of PCL clusters in gmt file and test set table do not match')

% ## Iterate for each PCL cluster to train OvR binary classifier using reference set similarity scores and PCL cluster-MOA agreement labels (1 for in-MOA / 0 for out-of-MOA) and apply to CGI profiles from unknown compounds to estimate their PCL cluster confidence scores

% save ss_tbl_reset as a copy of ss_tbl_sele which remains unchanged
% it will be used to subset and set ss_tbl_sele to the CGI profile similarity scores and MOA agreement labels for each specific PCL cluster iteratively 
ss_tbl_reset = ss_tbl_sele;

tic

out_tbl = {};
opt_tbl = {};
out_tbl_test_cmpd = {};

high_confidence_ppv_thres = high_confidence_pcl_confidence_score_thres; % defined as such for uniformity and simplicity across all PCL clusters

num_pcls = numel(pcl_list)

for i = 1:num_pcls

    loop_progress(i,num_pcls,100)

    choose_pcl = pcl_list(i);
    %choose_pcl
    
    jj = find(ismember(pcls_tbl.head, choose_pcl));
    %jj

    %pcls(jj).entry
    
    this_pcl_size = pcls(jj).len;
    %this_pcl_size

    % Subset to reference set rows for given PCL cluster %
    idx = strcmp(ss_tbl_reset.rid, choose_pcl);
    %sum(idx)

    ss_tbl_sele = ss_tbl_reset(idx,:);
    %size(ss_tbl_sele)

    choose_pcl_desc = ss_tbl_sele.pcl_desc(1);

    ss_tbl_sele = sortrows(ss_tbl_sele,{'pcl_similarity_score'},{'descend'});
    [unique_brds, idx] = unique(ss_tbl_sele.pert_id, 'stable'); % can be changed to broad_id, proj_broad_id, cid, etc. to consider top similarity scores at every compound detail level
    %numel(idx)

    % Commented out line would set highest_cmpd_min_pcl_similarity_score to similarity score of the reference set compound/condition with the lowest similarity to the PCL cluster to save compute time and storage
    % instead use -Inf to have confidence scores for similarity scores from all CGI profiles in the results file
    %highest_cmpd_min_pcl_similarity_score = min(ss_tbl_sele.pcl_similarity_score(idx))
    highest_cmpd_min_pcl_similarity_score = -Inf;

    %%% Take similarity scores to this PCL cluster for all reference set CGI profiles %%%
    
    all_trt_score_tbl = ss_tbl_sele(:, {'cid','pert_id','broad_id','proj_broad_id','target_description','pcl_similarity_score','rcorr','pcl_and_moa_agree'});

    all_trt_scores = all_trt_score_tbl.pcl_similarity_score(all_trt_score_tbl.pcl_similarity_score >= highest_cmpd_min_pcl_similarity_score);
    all_trt_labels = all_trt_score_tbl.pcl_and_moa_agree(all_trt_score_tbl.pcl_similarity_score >= highest_cmpd_min_pcl_similarity_score);

    % Take the highest similarity scoring conditions/CGI profiles for each reference set compound for initial binary classifier training to only count
    % compounds and not conditions as independent successes/true positives
    ss_tbl_sele = ss_tbl_sele(idx,:);
    
    %%% Section below uncommented would instead consider only highest the similarity scoring conditions for each reference set compound in the model fit %%% 
    %all_trt_scores = ss_tbl_sele.pcl_similarity_score;
    %all_trt_labels = ss_tbl_sele.pcl_and_moa_agree;
    
    %all_trt_score_tbl = ss_tbl_sele(:, {'cid','broad_id','target_description','pcl_similarity_score','rcorr','pcl_and_moa_agree'});
    
    %%%

    num_correct = sum(ss_tbl_sele.pcl_and_moa_agree == 1);
    %num_correct

    pcl_n = size(ss_tbl_sele, 1);

    tmp_opt_tbl = {};

    tmp_opt_tbl.pcl_rid = choose_pcl;
    tmp_opt_tbl.pcl_desc = choose_pcl_desc;
    tmp_opt_tbl.pcl_n = pcl_n;
    tmp_opt_tbl.pcl_n_correct = num_correct;
    tmp_opt_tbl.pcl_n_trt = length(all_trt_scores);
    tmp_opt_tbl.pcl_n_correct_trt = sum(all_trt_labels);

    if num_correct > 0

        train_cmpds = ss_tbl_sele.pert_id;
        
        train_pert_ids = ss_tbl_sele.pert_id;

        full_scores = ss_tbl_sele.pcl_similarity_score;

        assert(issorted(all_trt_scores, 'descend')) % high to low PCL cluster similarity score 

        assert(issorted(full_scores, 'descend')) % high to low PCL cluster similarity score

        full_train_scores = full_scores;

        full_labels = ss_tbl_sele.pcl_and_moa_agree; % binary labels for if reference set compound has shared MOA with PCL cluster (1 for in-MOA / 0 for out-of-MOA)

        % ROC analysis / fitting of binary classifier for current PCL cluster to threshold similarity scores to it (based on reference set) and assign mechanism 
        % if an unknown compound's similarity score exceeds the threshold
        % binarized class labels (in-MOA vs out-of-MOA), PCL cluster similarity scores, Positive class is 1 (in-MOA), Negative class is 0 (out-of-MOA)
        % prior probabilities are empirical (default) based on in-MOA/out-of-MOA frequencies of PCL cluster's MOA in reference set
        % default [0 1;1 0] (since MATLAB R2022a) misclassification costs specified
        [fpr,tpr,t,auc,optrocpt] = perfcurve(full_labels,full_train_scores,1,'NegClass',0,'UseNearest','off','XVals','all','Cost',[0 1; 1 0]);
        %auc
        %optrocpt
        % PPV (Positive Predictive Value) taken as empirical PCL cluster confidence scores at each observed PCL cluster similarity score from reference set
        % full_train_scores contains the highest similarity scoring CGI profile each reference set compound (pert_id)
        [ppv,tpr,t] = perfcurve(full_labels,full_train_scores,1,'XCrit','ppv','NegClass',0,'UseNearest','off','XVals','all','Cost',[0 1; 1 0]);

        % Extend PCL cluster similarity score to confidence score mapping to all reference set CGI profiles contained in all_trt_scores
        % through next-neighbor interpolation which assigns the PCL cluster confidence score of an unobserved (not maximal) CGI profile for a reference set compound
        % to the PCL cluster confidence score of the next most similar observed (maximal for the same or different reference set compound) CGI profile
        ppv_full = interp1(unique(full_scores, 'stable'), ppv(2:end), all_trt_scores,'next');
        
        % Ensure no missing values
        ppv_full = fillmissing(ppv_full, 'nearest');

        % Optimal threshold based on ROC analysis and default misclassification cost
        isOptimalPoint = ismember([fpr, tpr], optrocpt, 'rows');

        % Use find to return the indices of the rows that match the 1-by-2 array
        optimalPoint = find(isOptimalPoint(2:end)) + 1;

        tmp_opt_tbl.pcl_auc = auc;
        
        % high_confidence_ppv_thres set to 1; defined as such for uniformity and simplicity across all PCL clusters
            
        full_model_high_confidence_score_thres = all_trt_scores(max(find(single(ppv_full) >= single(high_confidence_ppv_thres))));
        %full_model_high_confidence_score_thres
        
        if isempty(full_model_high_confidence_score_thres) % should never be the case given definition of PCL clusters as having some region with PCL cluster confidence score = 1 but for completeness and error-checking
            full_model_high_confidence_score_thres = nan;
            
            full_model_next_score_thres = all_trt_score_tbl.pcl_similarity_score(1);
            full_model_next_score_ppv = 0;
            
            full_model_next_score_idx = 1;  
        else
            full_model_next_score_thres = all_trt_scores(max(find(single(ppv_full) >= single(high_confidence_ppv_thres)))+1);
            full_model_next_score_ppv = ppv_full(max(find(single(ppv_full) >= single(high_confidence_ppv_thres)))+1);
            
            full_model_next_score_idx = max(find(single(ppv_full) >= single(high_confidence_ppv_thres)))+1;
        end
        
        tmp_opt_tbl.pcl_high_confidence_similarity_score_threshold = full_model_high_confidence_score_thres;
        tmp_opt_tbl.pcl_high_confidence_ppv_thres = high_confidence_ppv_thres;

        if length(optimalPoint)>0

            optimalPPV = ppv(optimalPoint);

            full_model_opt_score_thres = all_trt_scores(max(find(single(ppv_full) >= single(optimalPPV))));
            %full_model_opt_score_thres
            
            %high_confidence_ppv_thres
            %optimalPPV

            % Uncommented would instead check for the similarity score, PPV, and MOA of the next compound after the optimal ROC point not high-confidence similarity score threshold
            %full_model_next_score_thres = all_trt_scores(max(find(single(ppv_full) >= single(optimalPPV)))+1)
            %full_model_next_score_ppv = ppv_full(max(find(single(ppv_full) >= single(optimalPPV)))+1)
            
            %full_model_next_score_idx = max(find(single(ppv_full) >= single(optimalPPV)))+1;

            tmp_opt_tbl.pcl_opt_similarity_score_threshold = full_model_opt_score_thres;
            tmp_opt_tbl.pcl_opt_ppv = ppv(optimalPoint);
            tmp_opt_tbl.pcl_opt_tpr = tpr(optimalPoint);
            tmp_opt_tbl.pcl_opt_fpr = fpr(optimalPoint);
        else

            optimalPPV = nan;

            full_model_opt_score_thres = nan;
            
            %high_confidence_ppv_thres
            %optimalPPV

            %full_model_next_score_thres = nan
            %full_model_next_score_ppv = nan
            
            %full_model_next_score_idx = 1;
            
            tmp_opt_tbl.pcl_opt_similarity_score_threshold = nan;
            tmp_opt_tbl.pcl_opt_ppv = nan;
            tmp_opt_tbl.pcl_opt_tpr = nan;
            tmp_opt_tbl.pcl_opt_fpr = nan;
        end

    else % should never be the case that there are no positive samples (reference set CGI profiles with matching MOA to the PCL cluster) for any PCL cluster but for completeness and error-checking
        ppv_full = repmat(nan,[size(all_trt_score_tbl, 1) 1]);
        
        full_model_high_confidence_score_thres = nan;
            
        full_model_next_score_thres = all_trt_score_tbl.pcl_similarity_score(1);
        full_model_next_score_ppv = 0;

        full_model_next_score_idx = 1;
        
        tmp_opt_tbl.pcl_high_confidence_similarity_score_threshold = nan;
        tmp_opt_tbl.pcl_high_confidence_ppv_thres = nan;
        
        optimalPPV = nan;
        
        full_model_opt_score_thres = nan;

        tmp_opt_tbl.pcl_auc = nan;

        tmp_opt_tbl.pcl_opt_similarity_score_threshold = nan;
        tmp_opt_tbl.pcl_opt_ppv = nan;
        tmp_opt_tbl.pcl_opt_tpr = nan;
        tmp_opt_tbl.pcl_opt_fpr = nan;

    end
    
    all_trt_score_tbl = all_trt_score_tbl(all_trt_score_tbl.pcl_similarity_score >= highest_cmpd_min_pcl_similarity_score, :);

    full_model_next_score_moa = all_trt_score_tbl.target_description(full_model_next_score_idx);

    tmp_opt_tbl.next_score_moa = full_model_next_score_moa;
    
    all_trt_score_tbl.pcl_confidence_score = ppv_full; % PPV (Positive Predictive Value) taken as empirical PCL cluster confidence scores at each observed PCL cluster similarity score from reference set
    
    all_trt_score_tbl.rid = repmat(choose_pcl,[size(all_trt_score_tbl, 1) 1]);
    all_trt_score_tbl.pcl_desc = repmat(choose_pcl_desc,[size(all_trt_score_tbl, 1) 1]);


    out_tbl = [out_tbl; all_trt_score_tbl];

    tmp_opt_tbl = struct2table(tmp_opt_tbl);
    
    tmp_opt_tbl.full_model_next_score_thres = full_model_next_score_thres;
    tmp_opt_tbl.full_model_next_score_ppv = full_model_next_score_ppv;

    opt_tbl = [opt_tbl; tmp_opt_tbl];
    
    
    % test compounds section % 
    
    % Subset to test compound set rows for given PCL cluster %
    test_idx = strcmp(ss_tbl_test_cmpd.rid, choose_pcl);
    %sum(test_idx)

    ss_tbl_sele_test_cmpd = ss_tbl_test_cmpd(test_idx,:);
    %size(ss_tbl_sele_test_cmpd)

    ss_tbl_sele_test_cmpd = sortrows(ss_tbl_sele_test_cmpd,{'pcl_similarity_score'},{'descend'});

    %%% Take similarity scores to this PCL cluster for all test compound CGI profiles %%%
    
    all_trt_score_tbl_test_cmpd = ss_tbl_sele_test_cmpd(:, {'cid','pert_id','broad_id','proj_broad_id','target_description','pcl_similarity_score','rcorr','pcl_and_moa_agree'});
    
    all_trt_scores_test_cmpd = all_trt_score_tbl_test_cmpd.pcl_similarity_score;
    all_trt_labels_test_cmpd = all_trt_score_tbl_test_cmpd.pcl_and_moa_agree;
    
    % all_trt_scores and ppv_full are from previous training step - similarity scores and confidence scores for all reference set CGI profiles to the PCL cluster
    assert(length(all_trt_scores) == length(ppv_full), 'Length of reference set PCL cluster similarity scores does not match length of PCL cluster confidence scores');
    [score_thresholds, idx_unique] = unique(all_trt_scores, 'stable');
    
    % Estimate PCL cluster confidence scores for each test compound CGI profile based on their PCL cluster confidence score
    % by fitting to mapping learned from all reference set CGI profiles through linear interpolation which considers the distance (or difference)
    % an unknown CGI profile’s similarity score is to the two reference set CGI profile similarity scores it is between when estimating the confidence score
    test_ppv = interp1(score_thresholds, ppv_full(idx_unique), all_trt_scores_test_cmpd, 'linear');
    
    % Ensure no missing values
    missing_indices = find(isnan(test_ppv));
    num_missing_indices = length(missing_indices);
    
    if num_missing_indices > 0
    
        filled_indices = find(~isnan(test_ppv));
        first_filled_index = filled_indices(1);
    
        missing_max_indices = missing_indices(missing_indices < filled_indices(1));
        num_missing_max_indices = length(missing_max_indices);
    
        missing_min_indices = missing_indices(missing_indices > filled_indices(1));
        num_missing_min_indices = length(missing_min_indices);
        
        if num_missing_max_indices > 0
            test_ppv(missing_max_indices) = max(ppv_full(idx_unique));
        end
        
        if num_missing_min_indices > 0
            test_ppv(missing_min_indices) = min(ppv_full(idx_unique));
        end
    
    end
    
    %test_ppv = fillmissing(test_ppv, 'nearest');

    assert(sum(isnan(test_ppv)) == 0)
    
    all_trt_score_tbl_test_cmpd.pcl_confidence_score = test_ppv; % Post-test probability / linear interpolation of fit PCL cluster Similarity Score to PPV curve from previous steps
    
    all_trt_score_tbl_test_cmpd.rid = repmat(choose_pcl,[size(all_trt_score_tbl_test_cmpd, 1) 1]);
    all_trt_score_tbl_test_cmpd.pcl_desc = repmat(choose_pcl_desc,[size(all_trt_score_tbl_test_cmpd, 1) 1]);


    out_tbl_test_cmpd = [out_tbl_test_cmpd; all_trt_score_tbl_test_cmpd];

    
    % End test compounds section %

end

toc

i

assert(i == num_pcls, sprintf('Not all PCL clusters processed, only %d PCL clusters', i))

headt(out_tbl)

headt(opt_tbl)

headt(out_tbl_test_cmpd)

% ## Save results in tabular/.txt form

% Save output tables as .txt files

wtable(out_tbl,fullfile(outdir,out_tbl_savename));

wtable(opt_tbl,fullfile(outdir,opt_tbl_savename));

wtable(out_tbl_test_cmpd,fullfile(outdir,out_tbl_test_cmpd_savename));

% Compress .txt files and delete to free up storage space (the compressed tab-separated .tsv/.txt files will be read in automatically in R)

gzip(fullfile(outdir,out_tbl_savename));

gzip(fullfile(outdir,opt_tbl_savename));

gzip(fullfile(outdir,out_tbl_test_cmpd_savename));

delete(fullfile(outdir,out_tbl_savename));

delete(fullfile(outdir,opt_tbl_savename));

delete(fullfile(outdir,out_tbl_test_cmpd_savename));

disp('Size and head of PCL cluster similarity x confidence score results table for training (reference) set reference compounds:')
size(out_tbl)
head(out_tbl)

disp('Size and head of per-PCL high-confidence similarity score threshold results table based on training/OvR binary classifier fitting using reference set compounds:')
size(opt_tbl)
head(opt_tbl)

disp('Size and head of PCL cluster similarity x confidence score results table for test compounds:')
size(out_tbl_test_cmpd)
head(out_tbl_test_cmpd)

% ## Combine reference set and unknown PCL cluster confidence score results tables and save as GCT

tic
% Combine out_tbl and out_tbl_test_cmpd keeping only the specified columns
combined_tbl_to_join = [out_tbl(:, {'cid', 'rid', 'pcl_confidence_score'}); out_tbl_test_cmpd(:, {'cid', 'rid', 'pcl_confidence_score'})];
size(combined_tbl_to_join)

% Perform an inner join with ss_tbl_train_and_test to maintain original row ordering
[out_tbl_train_and_test_cmpd, ileft] = innerjoin(ss_tbl_train_and_test, combined_tbl_to_join, 'Keys', {'cid', 'rid'});

[ileft_sorted, idx_to_sort_out_tbl_train_and_test_cmpd] = sort(ileft);

out_tbl_train_and_test_cmpd = out_tbl_train_and_test_cmpd(idx_to_sort_out_tbl_train_and_test_cmpd, :);

toc

disp('Size and head of PCL cluster similarity x confidence score results table for training (reference) and test set reference compounds:')
size(out_tbl_train_and_test_cmpd)
headt(out_tbl_train_and_test_cmpd)

% ss_tbl_train_and_test.cid(1:10)

% out_tbl_train_and_test_cmpd.cid(1:10)

% ss_tbl_train_and_test.rid(1:10)

% out_tbl_train_and_test_cmpd.rid(1:10)

disp('Length of ss_tbl_train_and_test cid column')
length(ss_tbl_train_and_test.cid)
disp('Length of out_tbl_train_and_test_cmpd cid column')
length(out_tbl_train_and_test_cmpd.cid)

assert(isequal(ss_tbl_train_and_test.cid, out_tbl_train_and_test_cmpd.cid), 'Ordering of CGI profiles (cids) in combined reference set and unknown compound results table differs from original ss_tbl_train_and_test table ordering')

assert(isequal(ss_tbl_train_and_test.rid, out_tbl_train_and_test_cmpd.rid), 'Ordering of PCL clusters (rids) in combined reference set and unknown compound results table differs from original ss_tbl_train_and_test table ordering')

% ### Section to prepare and save GCT

out_tbl_train_and_test_cmpd.pcl_id = out_tbl_train_and_test_cmpd.rid;

disp('Size and head of out_tbl_train_and_test_cmpd tbl:');
disp(size(out_tbl_train_and_test_cmpd));
disp(headt(out_tbl_train_and_test_cmpd));

% Make gctx
out_tbl_train_and_test_cmpd.cid_idx = grp2idx(out_tbl_train_and_test_cmpd.cid);
out_tbl_train_and_test_cmpd.pcl_id_idx = grp2idx(out_tbl_train_and_test_cmpd.pcl_id);

% Find indices
[~,ridx] = unique(out_tbl_train_and_test_cmpd.pcl_id_idx);
rid = out_tbl_train_and_test_cmpd.pcl_id(ridx);
[~,cidx] = unique(out_tbl_train_and_test_cmpd.cid_idx);
cid = out_tbl_train_and_test_cmpd.cid(cidx);

% Create metadata

col_names = out_tbl_train_and_test_cmpd.Properties.VariableNames;

% Find all column names that do not contain "pcl"
non_pcl_col_names = col_names(~contains(col_names, 'pcl') & ~contains(col_names, {'rid','cid_idx'}));
pcl_col_names = col_names(contains(col_names, 'pcl') & ~contains(col_names, {'pcl_similarity_score','pcl_and_moa_agree','pcl_id_idx','rank_pcl','pcl_confidence_score'}));

col_meta = out_tbl_train_and_test_cmpd(cidx, non_pcl_col_names);

disp('Size and head of final col_meta for gct files:');
disp(size(col_meta));
disp(headt(col_meta));

row_meta = out_tbl_train_and_test_cmpd(ridx, pcl_col_names);

disp('Size and head of final row_meta for gct files:');
disp(size(row_meta));
disp(headt(row_meta));

% Create an empty matrix
mat = nan(numel(rid),numel(cid));
ind = sub2ind(size(mat),out_tbl_train_and_test_cmpd.pcl_id_idx,out_tbl_train_and_test_cmpd.cid_idx);

ds = mkgctstruct(mat,'rid',rid','cid',cid);
ds = annotate_ds(ds,table2struct(col_meta),'dim','column','keyfield','cid');
ds = annotate_ds(ds,table2struct(row_meta),'dim','row','keyfield','pcl_id');

% Create gcts with various data
fields = {'pcl_confidence_score'};

for ii = 1:numel(fields)
    disp(fields{ii})
    ds_tmp = ds;
    ds_tmp.mat(ind) = out_tbl_train_and_test_cmpd.(fields{ii});

    % Verify
    rid1 = ds_tmp.rid(2);
    cid1 = ds_tmp.cid(1);
    corr1 = ds_tmp.mat(2,1);
    idx = ismember(out_tbl_train_and_test_cmpd.cid,cid1)&ismember(out_tbl_train_and_test_cmpd.pcl_id,rid1);
    assert(sum(idx)>0,'Rows or columns of the gct structure do not match');
    assert(out_tbl_train_and_test_cmpd.(fields{ii})(idx)==corr1,'Value in the matrix do not match the one in the table');
    
    %ds_tmp = ds_slice(ds_tmp, 'rid', s.rid, 'cid', s.cid);
    
    % Ensure order of CGI profiles (cid) and PCL clusters (rid) are the same as input PCL cluster similarity score gct
    assert(isequal(ds_tmp.cid, s.cid),'Columns are not sorted identically')
    assert(isequal(ds_tmp.rid, s.rid(1:i)),'Rows are not sorted identically') % i used here from for-loop iterating over PCL clusters for testing purposes (if smaller i than num_pcls was used)
    
    mkgctx(fullfile(outdir,[fields{ii},'.gctx']), ds_tmp)
    %mkgct(fullfile(outdir,[fields{ii},'.gct']), ds_tmp, 'precision', 8)
    
end