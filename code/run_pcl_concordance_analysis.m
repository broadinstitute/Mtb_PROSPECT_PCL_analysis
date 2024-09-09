function run_pcl_concordance_analysis(wkdir,c_corr,c_corr_rank,gmt,makefig,figdir,makegct,gctdir)

mk_cd_dir(wkdir,true)
mk_cd_dir(figdir,false)
mk_cd_dir(gctdir,false)

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
	out(ii).geneset_size = gmt(ii).len;

	ridx = ismember(c_corr.rid, gmt(ii).entry);
	cidx = ismember(c_corr.cid, gmt(ii).entry);
	if sum(ridx)>0 & sum(cidx)>0
		try
			c_corr_tmp = hclust(ds_slice(c_corr,'ridx',ridx,'cidx',cidx),'is_pairwise',true);
			c_corr_rank_tmp = ds_slice(c_corr_rank,'rid',c_corr_tmp.rid,'cid',c_corr_tmp.cid);
		catch ME
			disp(ME)
			c_corr_tmp = ds_slice(c_corr,'ridx',ridx,'cidx',cidx);
			c_corr_rank_tmp = ds_slice(c_corr_rank,'rid',c_corr_tmp.rid,'cid',c_corr_tmp.cid);
		end

		if makegct
			mkgct(fullfile(wkdir,gctdir,[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr.gct']),c_corr_tmp)
			mkgct(fullfile(wkdir,gctdir,[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr_rank.gct']),c_corr_rank_tmp)
		end

		out(ii).median_corr = median_of_medians(c_corr_tmp.mat,'row',true);
		out(ii).median_corr_rank = median_of_medians(c_corr_rank_tmp.mat,'row',true);
		if makefig
			figure
				set(gcf,'PaperPosition',[0,0,30,26])
				imagesc(c_corr_tmp.mat)
				colormap(jet)
				caxis([-1,1])
				colorbar
				set(gca,'FontSize',14)
				set(gca,'XTick',1:sum(ridx))
				set(gca,'YTick',1:sum(ridx))
				set(gca,'TickLabelInterpreter','none')
				set(gca,'YTickLabel',c_corr.rid(ridx))
				set(gca,'XTickLabel',[])
				set(gca,'FontSize',12)
				title(sprintf('Pearson correlation for: %s <%.1f>, <%.1f>', out(ii).head{1}, out(ii).median_corr, out(ii).median_corr_rank),'FontSize',30)
			saveas_png(gcf,fullfile(wkdir,figdir),[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr.png'])
			close(gcf)
			
			figure
				set(gcf,'PaperPosition',[0,0,30,26])
				imagesc(c_corr_rank_tmp.mat)
				colormap(jet)
				caxis([1,160])
				colorbar
				set(gca,'FontSize',14)
				set(gca,'XTick',1:sum(ridx))
				set(gca,'YTick',1:sum(ridx))
				set(gca,'TickLabelInterpreter','none')
				set(gca,'YTickLabel',c_corr.rid(ridx))
				set(gca,'XTickLabel',[])
				set(gca,'FontSize',12)
				title(sprintf('Pearson correlation rank for: %s <%.1f>, <%.1f>', out(ii).head{1}, out(ii).median_corr, out(ii).median_corr_rank),'FontSize',30)
			saveas_png(gcf,fullfile(wkdir,figdir),[strrep(strrep(strrep(strrep(out(ii).head{1},' ','_'),'|','_'),'/','-'),':','-'),'_corr_rank.png'])
			close(gcf)
		end
	else
		out(ii).median_corr = nan;
		out(ii).median_corr_rank = nan;
	end
end

disp('Save output table')
out = struct2table(out)
wtable(out,fullfile(wkdir,[figdir,'_corr_summary.txt']))


