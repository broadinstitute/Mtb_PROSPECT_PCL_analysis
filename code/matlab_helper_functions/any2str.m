function out = any2str(v)

% Convert an array of numbers or a cell array of numbers 
% and strings to a cell array of strings  

out = v;
if ismember(class(out),{'single','double','logical'})==1
    out = arrayfun(@(s) num2str(s), out, 'uni',0);
elseif ismember(class(out),{'char','cell'})==1
    idx = ismember(cellfun(@class,out,'uni',0),{'double','logical'});
    if sum(idx)>0
        out(idx) = cellfun(@(s) num2str(s), out(idx), 'uni',0);
    end
else
    true;
end

