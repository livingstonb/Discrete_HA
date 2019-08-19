classdef Prefs < handle

	properties (SetAccess = private)
		betadist;
		betatrans;
		betagrid0;

	end

	methods
		function obj = initialize_discount_factor(params)
			% discount factor distribution
		    if  params.nb == 1
		        obj.betadist = 1;
		        obj.betatrans = 1;
		    elseif params.nb > 1
		        % Equal probability in stationary distribution
		        obj.betadist = ones(params.nb,1) / params.nb; 
		        % Probability of switching from beta_i to beta_j, for i=/=j
		        betaswitch_ij = params.betaswitch / (params.nb-1);
		        % Create matrix with (1-betaswitch) on diag and betaswitch_ij
		        % elsewhere
		        diagonal = (1-params.betaswitch) * ones(params.nb,1);
		        off_diag = betaswitch_ij * ones(params.nb);
		        off_diag = off_diag - diag(diag(off_diag));
		        obj.betatrans = off_diag + diag(diagonal);
		    end
		    obj.betacumdist = cumsum(obj.betadist);
		    obj.betacumtrans = cumsum(obj.betatrans,2);
		    
		    % Create grid - add beta to grid later since we may iterate
		    bw = params.betawidth;
		    switch p.nb
		        case 1
		            obj.betagrid0 = 0;
		        case 2
		            obj.betagrid0 = [-bw/2 bw/2]';
		        case 3
		            obj.betagrid0 = [-bw 0 bw]';
		        case 4
		            obj.betagrid0 = [-3*bw/2 -bw/2 bw/2 3*bw/2]';
		        case 5
		            obj.betagrid0 = [-2*bw -bw 0 bw 2*bw]';
		    end

	end

end