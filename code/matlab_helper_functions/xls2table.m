
function out = xls2table(inp, sheet, hasheader)
% out = xls2table(inp, sheet, hasheader)
% Function to parse an Excel spreadsheet and convert it into a table.

if ~ischar(inp)
    error('File name should be char')
end

assert(islogical(hasheader), 'Argument hasheader should be logical')

if exist(inp, 'file')
    [~,~,out] = xlsread(inp,sheet);
    if hasheader
        out(:,cellfun(@ischar, out(1,:))==0) = [];

        out = cell2table(out(2:end,:),'VariableNames',...
            matlab.lang.makeValidName(out(1,:)));
    else
        out = cell2table(out);
    end
    varnames = out.Properties.VariableNames;
    for ii = 1:numel(varnames)
        disp([num2str(ii), ' ', varnames{ii}])

        if iscell(out.(varnames{ii}))
            idx_nan = strcmp('NaN',out.(varnames{ii})) | ...
                strcmp('nan',out.(varnames{ii})) | ...
                strcmp('#N/A',out.(varnames{ii})) | ...
                strcmp('NaN',out.(varnames{ii}));
            out.(varnames{ii})(idx_nan) = {nan};
            
            idx_numeric = cellfun(@isnumeric, (out.(varnames{ii})));
            idx_char = cellfun(@ischar, (out.(varnames{ii})));
            
            if sum(idx_numeric)>0 && sum(idx_char)>0
                out.(varnames{ii})(idx_numeric) = ...
                    cellfun(@num2str, out.(varnames{ii})(idx_numeric),'uni',0);
            elseif sum(idx_numeric)==size(out,1)
                    out.(varnames{ii}) = cell2mat(out.(varnames{ii}));
            end
        end
    end
else
    error('File %s could not be found',inp) 
end
