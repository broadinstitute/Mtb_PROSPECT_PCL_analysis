function run_concordance_analysis(wkdir,c_corr,c_corr_rank,gmt,makefig,figdir,makegct,gctdir,show_hclust)

% Function to subset, save as individual gctx files, 
% and plot the correlation and rank of correlation matrices 
% for each MOA or treatment cluster (from spectral clustering)
% saves a '*_corr_summary.txt' file 
% 
% by default, average rank in correlation (symmetrized between pairs of treatments) is highlighted in the range of 1 to 160
%
% c_corr - path to or correlation matrix gct object
% c_corr_rank - path to or rank of correlation matrix gct object
% gmt - path to or gmt object with MOA or treatment cluster information
% makefig - boolean to make figures
% figdir - directory to save figures
% makegct - boolean to save gct files (in addition to gctx)
% gctdir - directory to save gctx/gct files
% show_hclust - boolean to apply and show hierarchical clustering of matrices in saved figures

mk_cd_dir(wkdir,false)
mk_cd_dir(fullfile(wkdir, figdir),false)
mk_cd_dir(fullfile(wkdir, gctdir),false)

if ischar(gmt)
	gmt = parse_gmt(gmt)
end

disp('Parse gctx with correlation')
if isstr(c_corr)
    c_corr = parse_gctx(c_corr);
end
% Make sure this file is symmetric
c_corr = ds_slice(c_corr,'rid',c_corr.cid);

disp('Parse gctx with correlation ranks')
if isstr(c_corr_rank)
    c_corr_rank = parse_gctx(c_corr_rank);
end
% Make sure this file is symmetric
c_corr_rank = ds_slice(c_corr_rank,'rid',c_corr.rid,'cid',c_corr.cid);

disp('Remove classes with less than 2 members')
gmt = gmt([gmt.len]>2);

disp('Run summary and make figures')
for ii = 1:numel(gmt)
	loop_progress(ii, numel(gmt), 20)
	out(ii).head = {strrep(gmt(ii).head,'"','')};
	out(ii).desc = {gmt(ii).desc};
	out(ii).size = gmt(ii).len; % number of dsCGI profiles/treatments

	ridx = ismember(c_corr.rid, gmt(ii).entry);
	cidx = ismember(c_corr.cid, gmt(ii).entry);
	if sum(ridx)>0 & sum(cidx)>0
		if show_hclust
			try
				c_corr_tmp = hclust(ds_slice(c_corr,'ridx',ridx,'cidx',cidx),'is_pairwise',true);
				c_corr_rank_tmp = ds_slice(c_corr_rank,'rid',c_corr_tmp.rid,'cid',c_corr_tmp.cid);
			catch ME
				disp(ME)
				c_corr_tmp = ds_slice(c_corr,'ridx',ridx,'cidx',cidx);
				c_corr_rank_tmp = ds_slice(c_corr_rank,'rid',c_corr_tmp.rid,'cid',c_corr_tmp.cid);
			end
		else
			c_corr_tmp = ds_slice(c_corr,'ridx',ridx,'cidx',cidx);
			c_corr_rank_tmp = ds_slice(c_corr_rank,'rid',c_corr_tmp.rid,'cid',c_corr_tmp.cid);
		end

		% default precision for saving gct files is 4 decimals
		if makegct
			mkgct(fullfile(wkdir,gctdir,[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr.gct']),c_corr_tmp); % default precision is 4 decimals
			mkgct(fullfile(wkdir,gctdir,[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr_rank.gct']),c_corr_rank_tmp); % default precision is 4 decimals
		end

		% save gctx files by default
		mkgctx(fullfile(wkdir,gctdir,[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr.gctx']),c_corr_tmp); % no loss in decimal precision
		mkgctx(fullfile(wkdir,gctdir,[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr_rank.gctx']),c_corr_rank_tmp); % no loss in decimal precision
        
		% Create blue-white-red colormap
		bwr_colormap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); ...
		ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];

		out(ii).median_corr = median_of_medians(c_corr_tmp.mat,'row',true);
		out(ii).median_corr_rank = median_of_medians(c_corr_rank_tmp.mat,'row',true);
		if makefig
			figure
				set(gcf,'PaperPosition',[0,0,30,26])
				imagesc(c_corr_tmp.mat)
				colormap(bwr_colormap)
				caxis([-1,1])
				colorbar
				set(gca,'FontSize',14)
				set(gca,'XTick',1:sum(ridx))
				set(gca,'YTick',1:sum(ridx))
				set(gca,'TickLabelInterpreter','none')
				set(gca,'YTickLabel',c_corr.rid(ridx))
				set(gca,'XTickLabel',[])
				set(gca,'FontSize',12)
				title(sprintf('Pearson correlation for: %s\nMedian correlation:<%.1f>\nMedian rank in correlation across KABX:<%.1f>', out(ii).head{1}, out(ii).median_corr, out(ii).median_corr_rank),'FontSize',30)
			saveas_png(gcf,fullfile(wkdir,figdir),[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr.png'])
			close(gcf)
			
			figure
				set(gcf,'PaperPosition',[0,0,30,26])
				imagesc(c_corr_rank_tmp.mat)
				colormap(hot)
				caxis([0,160])
				colorbar
				set(gca,'FontSize',14)
				set(gca,'XTick',1:sum(ridx))
				set(gca,'YTick',1:sum(ridx))
				set(gca,'TickLabelInterpreter','none')
				set(gca,'YTickLabel',c_corr.rid(ridx))
				set(gca,'XTickLabel',[])
				set(gca,'FontSize',12)
				title(sprintf('Average rank in Pearson correlation across KABX: %s\nMedian correlation:<%.1f>\nMedian rank in correlation across KABX:<%.1f>', out(ii).head{1}, out(ii).median_corr, out(ii).median_corr_rank),'FontSize',30)
			saveas_png(gcf,fullfile(wkdir,figdir),[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr_rank.png'])
			close(gcf)
		end
	else
		out(ii).median_corr = nan;
		out(ii).median_corr_rank = nan;
	end
end

disp('Save output table')
out = struct2table(out);
wtable(out,fullfile(wkdir,[figdir,'_corr_summary.txt']));

disp('Showing first 20 rows of the output table')
head(out, 20)

end