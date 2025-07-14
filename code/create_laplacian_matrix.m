function [Ln,k,en,den,k_gap_den,k_med_gap_den,k_num_zero_plus_one,k_num_zero] = create_laplacian_matrix(A,norm_method,make_fig,k_type)

	% Function to create a graph Laplacian L = D - A from an adjacency matrix A.
	% Laplacian can be used as is (norm_method = "none") or 
	% normalized using Ng-Jordan-Weiss method (norm_medthod = "symmetric").
	% K can be selected using several different approaches/eigengap heuristics:
	% - the largest gap in eigenvalues
	% - the index of the smallest non-zero eigenvalue
	% - the median (average) of the two above
	% - the number of eigenvalues equal to zero (or approximately, currently arbitrarily set to 100*eps)
	% [LN,K,EN,DEN,K_GAP_DEN,K_MED_GAP_DEN,K_NUM_ZERO_PLUS_ONE,K_NUM_ZERO] = CREATE_LAPLACIAN_MATRIX(A,NORM_METHOD,MAKE_FIG,K_TYPE)
    
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
den_mask(en >= 1) = 0; % exclude eigenvalues greater than or equal to one which would maximize rather than minimize the variance explained by clusters or connected components in the graph Laplacian derived from the MOA adjacency matrix; otherwise only all singleton CGI profile clusters could sometimes be selected for k_gap_den
[m,k_gap_den] = max(den_mask);

% Find the number of small eigenvalues approximately equal to zero
k_num_zero = sum(en<100*eps);

% Find the index of the smallest non-zero eigenvalue
k_num_zero_plus_one = k_num_zero + 1;

% Find average of the index of the smallest non-zero eigenvalue and the eigenvalue where the gap or change in consecutive eigenvalues was largest
k_med_gap_den = ceil(median([k_num_zero_plus_one k_gap_den])); % ensure whole number of clusters K selected using ceiling

% Select which of the approaches for estimating K to use
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
		error('Unknown k_type; Only "k_gap_den", "k_med_gap_den", "k_num_zero_plus_one", or "k_num_zero" can be used')
	end
end

if make_fig
	figure
	    plot(en,'o-', 'DisplayName','Eigenvalue')
	    hold on 
	    plot(den,'o-', 'DisplayName',sprintf('Difference in\nconsecutive\neigenvalues'))
	    plot([k,k],ylim, 'DisplayName',sprintf('Selected k'))
		yl = ylim;
		text(k+0.5,yl(2)*0.95,sprintf('k: %d',k),'interpreter','none','FontSize',12)
		legend('show', 'Location', 'east', 'FontSize', 12, 'TextColor', 'black', 'Box', 'off')
		xlabel(sprintf('Number of eigenvectors\n(sorted from lowest to highest eigenvalue)'))
end