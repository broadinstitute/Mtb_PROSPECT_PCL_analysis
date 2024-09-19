function out = find_match(s,e,is_match)

% OUT = FIND_MATCH(S,E,IS_MATCH)
% Function that looks for expression (E) in a string (S).
% If is_match is set to true, function returns 1 if E matches S
% and 0 if it doesn't.
% If is_match is set to false, function returns 1 if E doesn't
% match S and 0 if it does.

%disp(['regexp(s,',stringify(e),')'])
idx = regexp(s,stringify(e));
if is_match
	out = cellfun(@(x) ~isempty(x), idx);
else
	out = cellfun(@(x) isempty(x), idx);
end
