function [out_nonsingle,out,gmt] = spectral_clustering_for_pcls(c,cr,moa_class,thrsh,k_type,outdir,save_out,save_fig,random_seed)
	
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
    
	% Set random_seed number if not provided in the input
	if isempty(random_seed)
		random_seed = 0;
	end
	rng(random_seed, 'twister'); % set seed and rng for reproducibility; rng(0, 'twister') is the Matlab factory default

	% Parse gctx with correlation matrix if not provided in the input
	if isstr(c)
		c = parse_gctx(c);
	end
	
	% Parse gctx with correlation rank matrix if not provided in the input
	if isstr(cr)
		cr = parse_gctx(cr);
	end
	
	% Set threshold for the rank if not provided in the input
	if isempty(thrsh)
		%thrsh = numel(cr.rid)*1;
		thrsh = 20;
	end
	
	% Create adjacency matrix
	adj_mat = double(cr.mat<=thrsh);

	% Normalized Laplacian
	[L,k] = create_laplacian_matrix(adj_mat,'symmetric',true,k_type);

	figure
	    set(gcf,'PaperPosition',[0,0,20,40])
	    sgtitle(moa_class,'FontWeight','bold','FontSize',20,'Interpreter','none')
    
	    s1 = subplot(5,2,1);
	    imagesc(c.mat)
	    colormap(s1,jet)
	    caxis([-1,1])
	    colorbar
	    title('Pearson correlation')
    
	    s2 = subplot(5,2,2);
	    imagesc(cr.mat)
	    colormap(s2,flip(parula))
	    %caxis([0,thrsh*2])
		caxis([0,thrsh+1])
	    colorbar
	    title('Average rank in Pearson correlation across KABX')

		s3 = subplot(5,2,3);
		imagesc(adj_mat)
		colormap(s3, [1,1,1; 1,0,0.5])  % 1 for connected treatments, 0 for unconnected treatments
		caxis([0, 1])
		colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'})  % Set colorbar ticks and labels
		title(sprintf('Adjacency matrix\n(thresholded at %d)', thrsh))
    
	    s4 = subplot(5,2,4);
	    imagesc(L)
	    colormap(s4,jet)
	    caxis([min(L,[],'all'),max(L,[],'all')])
	    colorbar
	    title('Laplacian matrix\n(normalized via Ng-Jordan-Weiss method)')

	[idx,eigenvec,eigenval] = spectralcluster(adj_mat,k,'Distance','precomputed','LaplacianNormalization','symmetric');
	tmp = repmat(idx,1,size(adj_mat,2));
	tmp_plot = tmp;
	tmp = tmp==tmp';
	tmp_g = graph(double(tmp),'omitselfloops');
    
	    s5 = subplot(5,2,5);
	    plot(tmp_g)
	    f = findobj(s5);
	    f(2).NodeFontSize = 6;
	    title(sprintf('Graphs from spectral clustering (k = %d)', k))
    
	    %s6 = subplot(5,2,6);
	    %imagesc(tmp)
	    %colormap(s6,[1,1,1; 1,0,0]) % 1 for treatments in cluster together, 0 for treatments in different clusters
	    %caxis([0,1])
	    %colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'})  % Set colorbar ticks and labels
	    %title(sprintf('Clusters from spectral clustering (k = %d)',k))

		s6 = subplot(5,2,6);
	    imagesc(tmp_plot)
	    colormap(s6, jet(k)) % different color for each cluster ID
	    caxis([1,k])
	    colorbar 
	    title(sprintf('Clusters from spectral clustering (k = %d)',k))
		
	% Create output tables
	col_meta_tmp = cell2table([cr.cid,cr.cdesc],'VariableNames',['cid';cr.chd]);

	cluster_id = conncomp(tmp_g)';
	out = table(cluster_id);

	%%% FIRST: Annotate the treatment/cid accurately before anything else
	out.cid = cr.rid;
	%%% Fix added by AB on 06/11/2024

	cluster_size = array2table(tabulate(cluster_id),'VariableNames',{'cluster_id','cluster_size','percent'});
	cluster_size.percent = [];
	out = outerjoin(out,cluster_size,'Keys','cluster_id','MergeKeys',true,'Type','left');
	out.moa_class = repmat({moa_class},size(out,1),1);
	out = join(out,col_meta_tmp,'Keys','cid');

	%%% AB add cluster membership error check on 06/11/2024 %
	%checkOutTable = sortrows(out(:, {'cid', 'cluster_id'}), {'cluster_id', 'cid'});

	%checkTable = table(cr.cid, cluster_id, 'VariableNames', {'cid', 'cluster_id'});

	%checkTable = sortrows(checkTable, {'cluster_id', 'cid'});

	%assert(isequal(checkOutTable, checkTable), 'The tables are not identical.');
	%%%

	idx = out.cluster_size>1;
	if sum(idx)>0
		out_nonsingle = out(idx,:);
		out_nonsingle.group_id = strcat(moa_class,':group',any2str(out_nonsingle.cluster_id));

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