function [out_nonsingle,out,gmt,Ln,k,en,den,k_gap_den,k_med_gap_den,k_num_zero_plus_one,k_num_zero] = spectral_clustering_for_pcls(c,cr,moa_class,thrsh,k_type,outdir,save_out,save_fig,random_seed,show_hclust)
	
	% Run spectral clustering on a correlation rank matrix for a given class of compounds
	% [OUT_NONSINGLE,OUT,GMT] = SPECTRAL_CLUSTERING_FOR_PCLS(C,CR,MOA_CLASS,THRSH,OUTDIR,SAVE_OUT,SAVE_FIG)
	%
	% Inputs:
	% c - correlation coefficients for all members of the MOA class (gct structure or string)
	% cr - pairwise average rank of correlation coefficients across KABX for all members of the MOA class (gct structure or string)
	% moa_class - MOA class (string)
	% thrsh - threshold to convert cr to adjacency matrix (numeric)
	% k_type - approach for estimating number of K clusters
	% outdir - output directory to save outputs (string)
	% save_out - variable if output tables should be saved (logical)
	% save_fig - variable if ouput figures should be saved (logical)
	% random_seed - specified seed for initializing the random number generator, Matlab factory default is the Mersenne Twister generator with seed 0; can be used to repeat spectral clustering for different random seeds (k-means++ algorithm involves random sampling/centroid initialization)
	% show_hclust - boolean to apply hierarchical clustering to correlation matrix prior to running spectral clustering, whether or not this is performed can minimally impact final cluster results especially for larger MOAs with greater number of K clusters due to stochasticity (randomized initialization) of k-means++ clustering; original results did apply hierarchical clustering but this step is not required for the method to be successful
	%
	% spectral clustering is performed using adjacency matrix created from the pairwise average rank of correlation coefficients
	% adjacency matrix is created by thresholding treatment pairs based on their average rank of correlation coefficients for thrsh (or less) mutual nearest-neighbors
	% since the adjacency matrix is used: 'Distance','precomputed' is set for spectral clustering
	% Ng-Jordan-Weiss method normalized Laplacian is computed from adjacency matrix: 'LaplacianNormalization','symmetric' is set for spectral clustering
	%
	% k-means++ clustering is performed for selected K once with random_seed set and default parameters:
	% n-by-k matrix of first K eigenvectors of the normalized Laplacian matrix
	% 'Display','off' (default)
	% 'Distance','sqeuclidean' Squared Euclidean distance (default)
	% 'MaxIter',100 (default)
	% 'EmptyAction','singleton' (default)
	% 'OnlinePhase','off' (default)

	% Parse gctx with correlation matrix if not provided in the input
	if isstr(c)
		c = parse_gctx(c);
	end
	
	% Parse gctx with correlation rank matrix if not provided in the input
	if isstr(cr)
		cr = parse_gctx(cr);
	end
	
    
	% Set show_hclust if not provided in the input
	if isempty(show_hclust)
		show_hclust = true; 
	end

	plot_spacer = 0;
	plot_row_spacer = 0;
	c_orig = c;
	cr_orig = cr;
	if show_hclust
		idx_target_description = find(strcmp(cr.chd, 'target_description'));
		idx_broad_id = find(strcmp(cr.chd, 'broad_id'));
		idx_pert_dose = find(strcmp(cr.chd, 'pert_dose'));
		brd_id_dose_tbl = cell2table([c.cid, c.cdesc(:, [idx_target_description, idx_broad_id, idx_pert_dose])], 'VariableNames', {'cid', 'target_description', 'broad_id', 'pert_dose'});
		brd_id_dose_tbl = sortrows(brd_id_dose_tbl,{'target_description','broad_id','pert_dose'});
		c = ds_slice(c,'rid',brd_id_dose_tbl.cid,'cid',brd_id_dose_tbl.cid); % sort correlation matrix by broad_id and ascending dose
		try
			c = hclust(c,'is_pairwise',true); % apply hierarchical clustering
		catch ME
			disp(ME)
		end
		cr = ds_slice(cr,'rid',c.rid,'cid',c.cid); % sort correlation rank matrix to match correlation matrix
		plot_spacer = 2;
		plot_row_spacer = 1;
	end
    
	% Set random_seed number if not provided in the input
	if isempty(random_seed)
		random_seed = 0; 
	end
	rng(random_seed, 'twister'); % set seed and rng for reproducibility; rng(0, 'twister') is the Matlab factory default
    
	% Set threshold for the rank if not provided in the input
	if isempty(thrsh)
		%thrsh = numel(cr.rid)*1;
		%thrsh = round(log(numel(cr.rid)) * 20);
		thrsh = 20;
	end
	
	% Create adjacency matrix
	adj_mat = double(cr.mat<=thrsh);

	% Normalized Laplacian and find eigenvalues and eigenvectors to estimate number of K clusters in MOA class
	[Ln,k,en,den,k_gap_den,k_med_gap_den,k_num_zero_plus_one,k_num_zero] = create_laplacian_matrix(adj_mat,'symmetric',true,k_type);
	%moa_class_str = strrep(strrep(moa_class,' ','_'),'/','-');
    %saveas_png(gcf,outdir,[moa_class_str,'_eigenvalues.png'])
    %close(gcf)

	% Create blue-white-red colormap
	bwr_colormap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); ...
	ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];

	figure
	    set(gcf,'PaperPosition',[0,0,20,40])
	    sgtitle(moa_class,'FontWeight','bold','FontSize',20,'Interpreter','none')
    
	    s1 = subplot(5+plot_row_spacer,2,1);
	    imagesc(c_orig.mat)
	    colormap(s1,bwr_colormap)
	    caxis([-1,1])
	    colorbar
	    title('Pearson correlation')
    
	    s2 = subplot(5+plot_row_spacer,2,2);
	    imagesc(cr_orig.mat)
		%colormap(s2,flip(hot))
	    colormap(s2,hot)
	    %caxis([1,thrsh*2])
		caxis([1,thrsh+1])
	    colorbar
	    title(sprintf('Average rank in Pearson\ncorrelation across KABX'))

		if show_hclust
		    s1_hc = subplot(5+plot_row_spacer,2,1+plot_spacer);
		    imagesc(c.mat)
		    colormap(s1_hc,bwr_colormap)
		    caxis([-1,1])
		    colorbar
		    title(sprintf('Pearson correlation\n(hierarchical clustering)'))
    
		    s2_hc = subplot(5+plot_row_spacer,2,2+plot_spacer);
		    imagesc(cr.mat)
			%colormap(s2_hc,flip(hot))
		    colormap(s2_hc,hot)
		    %caxis([1,thrsh*2])
			caxis([1,thrsh+1])
		    colorbar
		    title(sprintf('Average rank in Pearson\ncorrelation across KABX\n(re-sorted by hierarchical clustering)'))
		end

		s3 = subplot(5+plot_row_spacer,2,3+plot_spacer);
		imagesc(adj_mat)
		colormap(s3, [1,1,1; 1,0,0.5])  % 1 for connected treatments, 0 for unconnected treatments
		caxis([0, 1])
		colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'})  % Set colorbar ticks and labels
		title(sprintf('Adjacency matrix\n(thresholded at %d)', thrsh))
    
	    s4 = subplot(5+plot_row_spacer,2,4+plot_spacer);
	    imagesc(Ln)
	    colormap(s4,jet)
	    caxis([min(Ln,[],'all'),max(Ln,[],'all')])
	    colorbar
	    title(sprintf('Laplacian matrix\n(normalized via Ng-Jordan-Weiss method)'))

		s5 = subplot(5+plot_row_spacer,2,5+plot_spacer);
	    plot(en,'o-', 'DisplayName','Eigenvalue')
	    hold on 
	    plot(den,'o-', 'DisplayName',sprintf('Difference in\nconsecutive\neigenvalues'))
	    plot([k,k],ylim, 'DisplayName',sprintf('Selected k'))
		yl = ylim;
		text(k+0.5,yl(2)*0.95,sprintf('k: %d',k),'interpreter','none','FontSize',12)
		legend('show', 'Location', 'east', 'FontSize', 12, 'TextColor', 'black', 'Box', 'off')
		xlabel(sprintf('Number of eigenvectors\n(sorted from lowest to highest eigenvalue)'))

	[idx,eigenvec,eigenval] = spectralcluster(adj_mat,k,'Distance','precomputed','LaplacianNormalization','symmetric');
	tmp = repmat(idx,1,size(adj_mat,2));
	tmp_plot = tmp;
	tmp = tmp==tmp'; % symmetric, block, logical matrix where 1 means the two treatments are in the same cluster, 0 otherwise
	tmp_g = graph(double(tmp),'omitselfloops');

	% Create output tables
	col_meta_tmp = cell2table([cr.cid,cr.cdesc],'VariableNames',['cid';cr.chd]);
	
	cluster_id = conncomp(tmp_g)'; % get cluster numbers (ids) for each treatment, note: these are different numerically than idx (from spectralcluster) but the same in terms of cluster membership; the actual number/id is arbitrary
	out = table(cluster_id);

	%%% Annotate the treatment/cid
	out.cid = cr.rid;

	cluster_size = array2table(tabulate(cluster_id),'VariableNames',{'cluster_id','cluster_size','percent'});
	cluster_size.percent = [];
	out = outerjoin(out,cluster_size,'Keys','cluster_id','MergeKeys',true,'Type','left'); % outerjoin operation automatically re-sorts the rows by the key (cluster_id)
	out.moa_class = repmat({moa_class},size(out,1),1);
	out = join(out,col_meta_tmp,'Keys','cid');

	%%% Cluster membership error check - check cluster ids are accurately annotated to treatment after outerjoin operation %%%
	checkOutTable = sortrows(out(:, {'cid', 'cluster_id'}), {'cluster_id', 'cid'});

	checkTable = table(cr.cid, cluster_id, 'VariableNames', {'cid', 'cluster_id'});

	checkTable = sortrows(checkTable, {'cluster_id', 'cid'});

	assert(isequal(checkOutTable, checkTable), 'The tables are not identical.');
	%%%
	
	% Create block matrix with clusters aligned along the diagonal
	tmp_cluster_id = repmat(cluster_id,1,size(adj_mat,2)); % repeat the cluster_id vector rowwise to create a matrix
	
	tmp_cluster_id(~tmp) = 0; % apply the logical, block matrix to set the treatments that are not in the same cluster to 0 and keep the treatments that are in the same cluster as the numeric cluster_id
	
	[sortedLabels, sortOrder] = sort(cluster_id);

	tmp_cluster_id_sorted_by_cluster = tmp_cluster_id(sortOrder, sortOrder); % sort the block matrix by cluster_id so that clusters are aligned along the diagonal
	
	c_mat_sorted_by_cluster = c.mat(sortOrder, sortOrder); % sort the correlation matrix by cluster_id so that clusters are aligned along the diagonal
	
	cr_mat_sorted_by_cluster = cr.mat(sortOrder, sortOrder); % sort the correlation rank matrix by cluster_id so that clusters are aligned along the diagonal

	[orig_check, orig_indices] = ismember(out.cid, c_orig.cid);

	assert(all(orig_check), 'Problem sorting cluster membership matrix by original treatment order')
    
	tmp_cluster_id_sorted_by_orig_trt_order = tmp_cluster_id(orig_indices, orig_indices); % sort the block matrix by the original treatment order of the correlation matrix
    
		s6 = subplot(5+plot_row_spacer,2,6+plot_spacer);
		imagesc(tmp_cluster_id_sorted_by_orig_trt_order)
		colormap(s6, colorcube(k+1)) % different color for each cluster ID
	    %colormap(s6, lines(k+1)) % different color for each cluster ID
		caxis([0,k])
		colorbar 
		title(sprintf('Clusters from spectral clustering (k = %d)\n(original treatment order)',k))

		s7 = subplot(5+plot_row_spacer,2,7+plot_spacer);
		plot(tmp_g)
		f = findobj(s7);
		f(2).NodeFontSize = 6;
		if show_hclust
			title(sprintf('Graphs from spectral clustering (k = %d)\n(row index is hierarchical clustering treatment order)', k))
		else
			title(sprintf('Graphs from spectral clustering (k = %d)\n(row index is original treatment order)', k))
		end
		
		s8 = subplot(5+plot_row_spacer,2,8+plot_spacer);
		imagesc(tmp_cluster_id_sorted_by_cluster)
		colormap(s8, colorcube(k+1)) % different color for each cluster ID
	    %colormap(s8, lines(k+1)) % different color for each cluster ID
		caxis([0,k])
		colorbar 
		title(sprintf('Clusters from spectral clustering (k = %d)\n(re-sorted by cluster)',k))
		
		s9 = subplot(5+plot_row_spacer,2,9+plot_spacer);
		imagesc(c_mat_sorted_by_cluster)
		colormap(s9,bwr_colormap)
		caxis([-1,1])
		colorbar
		title(sprintf('Pearson correlation\n(re-sorted by cluster)'))
		
		s10 = subplot(5+plot_row_spacer,2,10+plot_spacer);
		imagesc(cr_mat_sorted_by_cluster)
		%colormap(s10,flip(hot))
		colormap(s10,hot)
		%caxis([1,thrsh*2])
		caxis([1,thrsh+1])
		colorbar
		title(sprintf('Average rank in Pearson\ncorrelation across KABX\n(re-sorted by cluster)'))


    out.group_id = strcat(moa_class,':group',any2str(out.cluster_id));
	idx = out.cluster_size>1; % only keep clusters with more than one treatment
	if sum(idx)>0
		out_nonsingle = out(idx,:);
		%out_nonsingle.group_id = strcat(moa_class,':group',any2str(out_nonsingle.cluster_id));

		% Create gmt file
		gmt = tbl2gmt(table2struct(out_nonsingle),'group_field','group_id','desc_field','moa_class','member_field','cid');
	else
		out_nonsingle = {};
		gmt = {};
	end
	
	moa_class_str = strrep(strrep(moa_class,' ','_'),'/','-');
	if save_fig
	    saveas_png(gcf,outdir,[moa_class_str,'_spectral_clust.png'])
	end

	if save_out
	    wtable(out,fullfile(outdir,[moa_class_str,'_spectral_clust.txt']))
	    wtable(out_nonsingle,fullfile(outdir,[moa_class_str,'_spectral_clust_nonsingle.txt']))
	    mkgmt(fullfile(outdir,[moa_class_str,'_spectral_clust.gmt']),gmt)
	end
end