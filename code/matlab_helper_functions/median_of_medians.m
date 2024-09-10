function [out,out_vector] = median_of_medians(mat,direc,nodiag)

% [OUT, OUT_VECTOR] = MEDIAN_OF_MEDIANS(MAT, DIREC, NODIAG)
% Function to calculate a median of medians along rows or columns
% MAT - matrix
% DIREC - 'row' or 'column'
% NODIAG - true or false whether to include diagonal or not

if nodiag
    if size(mat,1)==size(mat,2)
        mat = set_diag(mat,nan);
    else
        error('Non-squared matrix was provided')
    end
end

out = nan;
if isequal(size(mat), [2,2])
    out = nanmedian([mat(1,2),mat(2,1)]);
else 
    switch direc
        case 'row'
            out_vector = nanmedian(mat,2);
            out = median(out_vector);
        case 'column'
            out_vector = nanmedian(mat,1);
            out = median(out_vector);
    end
end