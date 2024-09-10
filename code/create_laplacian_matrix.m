function [Ln,k,en,den,k_gap_den,k_num_zero] = create_laplacian_matrix(A,norm_method,make_fig,k_type)

	% Function to create a laplacian L = D - A from an adjacency matrix A.
	% Laplacian can be used as is (norm_method = "none") or 
	% normalized using Ng-Jordan-Weiss method (norm_medthod = "symmetric").
	% K is selected using two different approaches:
	% - the gap in eigenvalues
	% - the number of eigenvalues with low values (currently arbitrarily set to 100*eps)
	% [LN,K,EN,DEN,K_GAP_DEN,K_NUM_ZERO] = CREATE_LAPLACIAN(A,NORM_METHOD,MAKE_FIG,K_TYPE)
    
D = zeros(size(A));
D = set_diag(D,sum(A,2));
L = D - A;
    
switch norm_method
	case 'none'
		Ln = L;
	case 'symmetric'
		Ln = D^-0.5*L*D^-0.5;
	otherwise
		error('Only allowed options for norm_method are: "none" or "symmetric"')
end
en = eig(Ln);

% Find the largest difference in eigenvalues
den = diff(en);
den_mask = den;
den_mask(en >= 1) = 0;
[m,k_gap_den] = max(den_mask);

% Find the number of small eigenvalues (very arbitrarily defined)
k_num_zero = sum(en<100*eps);
k_num_zero_plus_one = k_num_zero + 1;

k_med_gap_den = ceil(median([k_num_zero_plus_one k_gap_den]));

% Select which of the two k_gaps to use
if sum(en<=eps)==size(A,1);
	k = size(A,1);
else
	switch k_type
	case 'k_gap_den'
		k = k_gap_den;
	case 'k_med_gap_den'
		k = k_med_gap_den;
	case 'k_num_zero'
		k = k_num_zero;
	case 'k_num_zero_plus_one'
		k = k_num_zero_plus_one;
	otherwise
		error('Unknown k_type; Only "k_gap_den" or "k_num_zero" can be used')
	end
end

if make_fig
	figure
	    plot(en,'o-')
	    hold on 
	    plot(den,'o-')
	    plot([k,k],ylim)
		yl = ylim;
		text(k+5,yl(2)*0.95,sprintf('k: %d',k),'interpreter','none','FontSize',12)
end