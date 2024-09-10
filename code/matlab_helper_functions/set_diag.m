function out_mat = set_diag(inp_mat, val)

% OUT_MAT = SET_DIAG(INP_MAT, VAL)
% Function to set diagonal in inp_mat to val

n = size(inp_mat,1);
inp_mat(1:n+1:n^2) = val;

out_mat = inp_mat;