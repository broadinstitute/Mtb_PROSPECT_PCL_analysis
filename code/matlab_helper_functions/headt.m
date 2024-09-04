function out = headt(tbl,varargin)

% OUT = HEADT(TBL, N)
% Function to print transposed header and first N rows from a table tbl.

if nargin>1
    n = varargin{1};
    assert(isnumeric(n),'Second argument of this function should be numeric')
else
    n = 1;
end

assert(istable(tbl), 'Provided argument does not represent a table')

out = cell2table(gcthd2idx([tbl.Properties.VariableNames; table2cell(tbl(1:n,:))]'),...
	'VariableNames',{'idx','field','value'});
