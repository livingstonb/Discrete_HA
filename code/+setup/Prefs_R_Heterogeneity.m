classdef Prefs_R_Heterogeneity < handle
	% This class stores values, distributions, and
	% transition matrices for heterogeneity in
	% beta, returns, or another variable (z).
	%
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties (SetAccess = private)
		betadist;
		betatrans;
		betagrid0;
		betacumdist;
		betacumtrans;

		zdist;
		ztrans;
		zcumdist;
		zcumtrans;

		rdist;
		rtrans;
		rcumdist;
		rcumtrans;

	end

	methods
		function obj = Prefs_R_Heterogeneity(params)
			obj.initialize_discount_factor(params);
			obj.initialize_IES_heterogeneity(params);
            obj.initialize_returns_heterogeneity(params);
		end

		%% -------------------------------------------------------
	    % Discount Factor Heterogeneity
	    % --------------------------------------------------------
		function obj = initialize_discount_factor(obj, params)
			% discount factor distribution
		    if  params.nbeta == 1
		        obj.betadist = 1;
		        obj.betatrans = 1;
		    elseif params.nbeta > 1
		        % Equal probability in stationary distribution
		        if numel(params.beta_dist) == 1
		        	obj.betadist = ones(params.nbeta,1) / params.nbeta;
		        elseif (numel(params.beta_dist)==params.nbeta) && (sum(params.beta_dist)==1)
		        	obj.betadist = params.beta_dist(:);
		        else
		        	error('Invalid distribution for betas')
		        end

		        if (numel(params.beta_dist) > 1) && (params.betaswitch>0)
					error('Model does not allow for setting both a probability of beta switching and the stationary distribution for beta')
				end
		        % Probability of switching from beta_i to beta_j, for i=/=j
		        betaswitch_ij = params.betaswitch / (params.nbeta-1);
		        % Create matrix with (1-betaswitch) on diag and betaswitch_ij
		        % elsewhere
		        diagonal = (1-params.betaswitch) * ones(params.nbeta,1);
		        off_diag = betaswitch_ij * ones(params.nbeta);
		        off_diag = off_diag - diag(diag(off_diag));
		        obj.betatrans = off_diag + diag(diagonal);
		    end
		    obj.betacumdist = cumsum(obj.betadist);
		    obj.betacumtrans = cumsum(obj.betatrans,2);
		    
		    % Create grid - add beta to grid later since we may iterate
		    bw = params.betawidth;
		    switch params.nbeta
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

		%% -------------------------------------------------------
	    % IES Heterogeneity (Epstein-Zin only)
	    % --------------------------------------------------------
	    function obj = initialize_IES_heterogeneity(obj, params)
		    if numel(params.risk_aver) > 1 || ((numel(params.invies) > 1)...
		    	&& (params.EpsteinZin == 1))

		    	obj.zdist = ones(params.nb,1) / params.nb;
		    	zswitch_ij = params.IESswitch / (params.nb-1);

		    	diagonal = (1-params.IESswitch) * ones(params.nb,1);
		    	off_diag = zswitch_ij * ones(params.nb);
		    	off_diag = off_diag - diag(diag(off_diag));
		    	obj.ztrans = off_diag + diag(diagonal);
		        
		        obj.zcumdist = cumsum(obj.zdist);
		        obj.zcumtrans = cumsum(obj.ztrans,2);
		    else
		        obj.zdist = 1;
		        obj.ztrans = 0;
		        obj.zcumdist = 1;
		        obj.zcumtrans = 0;
		    end
		end

		%% -------------------------------------------------------
	    % Returns Heterogeneity
	    % --------------------------------------------------------
		function initialize_returns_heterogeneity(obj, params)
            nr = numel(params.r);
			if nr > 1
		        obj.rdist = ones(nr,1) / nr;
		        rswitch = 0;

		        diagonal = (1-rswitch) * ones(nr,1);
		        off_diag = rswitch * ones(nr);
		        off_diag = off_diag - diag(diag(off_diag));
		        obj.rtrans = off_diag + diag(diagonal);
		        
		        obj.rcumdist = cumsum(obj.rdist);
		        obj.rcumtrans = cumsum(obj.rtrans,2);
		    else
		        obj.rdist = 1;
		        obj.rtrans = 0;
		        obj.rcumdist = 1;
		        obj.rcumtrans = 0;
		    end
		end
	end

end