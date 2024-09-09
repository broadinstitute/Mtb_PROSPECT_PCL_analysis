function out = gcthd2idx(hd)
% Function to print elements of column or row headers of a gct structure
% together with their indices
%
% m = gcthd2idx(ds.chd)

out = [num2cell(([1:length(hd)]')),hd];
