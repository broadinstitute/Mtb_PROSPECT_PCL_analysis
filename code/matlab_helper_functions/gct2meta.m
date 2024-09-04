function [row_meta,col_meta] = gctmeta2tbl(ds)

% [ROW_META, COL_META] = GCT2TBL(Ds)
% Extract gct to row and column metadata
% ds - input gct

% Load the table if provided as a string
if isstr(ds)
	disp(sprintf('%s> Loading input gct\n',mfilename))
	ds = parse_gctx(ds);
end

disp(sprintf('%s> Creating row and column metadata tables\n',mfilename))
row_meta = cell2table([ds.rid, ds.rdesc], 'VariableNames', ['rid'; ds.rhd]);
col_meta = cell2table([ds.cid, ds.cdesc], 'VariableNames', ['cid'; ds.chd]);
