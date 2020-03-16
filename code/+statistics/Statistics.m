classdef Statistics < handle

	properties
		pmf;
		pmf_a;
		cdf_a;

		beta_A;
		beta_Q;

		mean_a;
		mean_s;
		mean_x;
		median_a;

		mean_gross_y_annual;
		std_log_gross_y_annual;
		std_log_net_y_annual;
		annual_inc_dollars;

		wpercentiles;

		w_top10share;
		w_top1share;
		wgini;

		mpcs;

		constrained;
	end

	properties (Access=protected)
		p;
		income;
		grdDST;
		model;

		nx;
		freq;
	end

	methods
		function obj = Statistics(p, income, grdDST, model)
			obj.p = p;
			obj.income = income;
			obj.grdDST = grdDST;
			obj.model = model;

			obj.nx = p.nx_DST;
			obj.freq = p.freq;

			obj.pmf = model.adist;
		end

		function compute_statistics(obj)
			% obj.add_params();
			obj.compute_intro_stats();
			obj.construct_distributions();
			obj.compute_percentiles();
			obj.compute_inequality();
			% obj.compute_constrained();
			% obj.compute_deposit_stats();
		end

		function add_mpcs(obj, mpcs_obj)
			obj.mpcs = mpcs_obj.mpcs;

			% if mpc_obj.options.liquid_mpc
			% 	asset_indicator = 0;
			% 	mpc_type = '';
			% else
			% 	asset_indicator = 2;
			% 	mpc_type = 'illiq';
			% end
			% sfill2 = @(x,y) sfill(x, y, asset_indicator);

			% empty_stat = sfill2([], []);

		% 	empty_mpc_struct = struct(...
		% 		'shock_normalized', empty_stat,...
		% 		'shock', empty_stat,...
		% 		'quarterly', empty_stat,...
		% 		'annual', empty_stat...
		% 	);

		% 	nshocks = numel(obj.p.mpc_shocks);
		% 	for ishock = 1:nshocks
		% 		shock = obj.p.mpc_shocks(ishock);
		% 		shock_label = obj.p.quantity2label(shock);
		% 		mpcs_stats(ishock) = empty_mpc_struct;

		% 		mpcs_stats(ishock).shock_normalized = sfill2(...
		% 			shock * 100, 'Shock size, (% of mean ann inc)');

		% 		mpcs_stats(ishock).shock = sfill2(shock_label,...
		% 			'Shock size');

		% 		tmp = 100 * mpc_obj.mpcs(ishock).quarterly(1);
		% 		label = sprintf(...
		% 			'Quarterly %s MPC (%%), out of %s',...
		% 			mpc_type, shock_label);
		% 		mpcs_stats(ishock).quarterly = sfill2(tmp, label);

		% 		tmp = 100 * mpc_obj.mpcs(ishock).annual;
		% 		label = sprintf(...
		% 			'Annual %s MPC (%%), out of %s',...
		% 			mpc_type, shock_label);
		% 		mpcs_stats(ishock).annual = sfill2(tmp, label);
		% 	end

		% 	if mpc_obj.options.liquid_mpc
		% 		obj.mpcs = mpcs_stats;

		% 		obj.mpcs_over_ss = cell(1, nshocks);
		% 		for ishock = 1:nshocks
		% 			obj.mpcs_over_ss{ishock} =  mpc_obj.mpcs(ishock).mpcs(:,1);
		% 		end
		% 	else
		% 		obj.illiquid_mpcs = mpcs_stats;
		% 	end
		% end
		end
	end

	methods (Access=protected)
		function compute_intro_stats(obj)
			obj.beta_A = sfill(obj.p.beta0 ^ obj.freq,...
				'beta (annualized)');
		    obj.beta_Q = sfill(obj.p.beta0 ^ (obj.freq/4),...
		    	'beta (quarterly)');
		    
		    tmp = obj.expectation(obj.grdDST.a.matrix);
		    obj.mean_a = sfill(tmp, 'Mean wealth');

		    xdist = obj.model.xdist(:);

		    mean_s = dot(obj.model.sav_x(:), xdist);
		    obj.mean_s = sfill(mean_s, 'Mean s');

		    mean_x = dot(obj.model.xvals(:), xdist);
		    obj.mean_x = sfill(mean_x, 'Mean x');

		    % Income
		    mean_y = dot(obj.model.y_x(:) * obj.freq, xdist);
		    obj.mean_gross_y_annual = sfill(mean_y,...
		    	'Mean gross annual income');

		    if obj.freq == 1
		    	demeaned2 = (obj.model.y_x(:) - mean_y) .^ 2;
			    stdev_y = dot(demeaned2, xdist);
			    obj.std_log_gross_y_annual = sfill(stdev_y,...
			    	'Stdev log gross annual income');

			    demeaned2 = (obj.model.nety_x(:) - mean_y) .^ 2;
			    stdev_y = dot(demeaned2, xdist);
			    obj.std_log_net_y_annual = sfill(stdev_y,...
			    	'Stdev log net annual income');
			else
				obj.std_log_gross_y_annual = sfill(NaN,...
			    	'Stdev log gross annual income');
				obj.std_log_net_y_annual = sfill(NaN,...
			    	'Stdev log net annual income');
			end

			dollars = sprintf('$%g', obj.p.annual_inc_dollars);
			obj.annual_inc_dollars = sfill(dollars,...
			    'Dollar value of mean gross ann inc (numeraire)');
		end

		function construct_distributions(obj)
			obj.pmf_a = sum(reshape(obj.pmf, obj.nx, []), 2);
		    obj.cdf_a = cumsum(obj.pmf_a);
		end

		function compute_percentiles(obj)
			w_pct = pct_interp(obj.grdDST.a.vec, obj.cdf_a);

			npct = numel(obj.p.percentiles);
		    obj.wpercentiles = cell(1, npct);
			for ip = 1:npct
				pct_at = obj.p.percentiles(ip);

				tmp_b = sprintf('Wealth, %gth pctile', pct_at);
				obj.wpercentiles{ip} = sfill(...
					w_pct(pct_at/100), tmp_b);
			end
			obj.median_a = sfill(w_pct(0.5), 'Median wealth');
		end

		function compute_inequality(obj)
			% import HACTLib.aux.interp_integral_alt
			% import HACTLib.aux.unique_sort
			% import HACTLib.aux.direct_gini

			% Top liquid wealth shares
			cum_share = cumsum(obj.grdDST.a.vec .* obj.pmf_a);
			cum_share = cum_share / obj.mean_a.value;
			wshare_interp = griddedInterpolant(obj.cdf_a,...
				cum_share, 'pchip', 'nearest');

			tmp = 1 - wshare_interp(0.9);
			obj.w_top10share = sfill(tmp, 'Wealth, top 10% share');

			tmp = 1 - wshare_interp(0.99);
			obj.w_top1share = sfill(tmp, 'Wealth, top 1% share');
			
			% Gini coefficient
			tmp = aux.direct_gini(obj.grdDST.a.vec, obj.pmf_a);
			obj.wgini = sfill(tmp, 'Gini coefficient, wealth');
		end

		function compute_constrained(obj)
			% Constrained by fraction of mean ann inc
			constrained_interp = constrained_interp(...
	        	obj.grdDST.a.vec, obj.cdf_a);

		    neps = numel(obj.p.epsilon);
		    obj.constrained = cell(1, neps);
		    for ip = 1:neps
				htm = obj.p.epsilon(ip);

				tmp = constrained_interp(htm);
				obj.constrained{ip} = sfill(tmp,...
					sprintf('a <= %g', htm));
			end

			% Wealth / (quarterly earnings) < epsilon
			a_over_inc = obj.grdDST.a.vec ./ ...
				(obj.income.netymat_broadcast * (p.freq / 4));
		    a_over_inc = repmat(a_over_inc, [1, 1, 1, obj.p.nb, 1]);
		    pmf_AY = obj.pmf(:) * shiftdim(obj.income.yTdist, -1);
		    sorted_mat = sortrows([a_over_inc(:), pmf_AY(:)]);

		    cdf_AY = cumsum(sorted_mat(:,2));
		    vals = sorted_mat(:,1);

		    ay_interp = constrained_interp(vals, cdf_AY);

			obj.a_lt_ysixth = sfill(...
				ay_interp(1/6), 'a_i <= y_i / 6 (biweekly earnings)');
			obj.a_lt_ytwelfth = sfill(...
				ay_interp(1/12), 'a_i <= y_i / 12 (weekly earnings)');
		end

		function out = expectation(obj, vals)
			out = dot(obj.pmf(:), vals(:));
		end
	end
end

function out = sfill(value, label)
	out = struct(...
		'value', value,...
		'label', label...
	);
end

function interp_out = pct_interp(values, cdf_x)
	[cdf_x_u, iu] = unique(cdf_x(:), 'first');
	values_u = values(iu);

	if numel(cdf_x_u) >= 2
		interp_out = griddedInterpolant(...
			cdf_x_u, values_u, 'pchip', 'nearest');
	else
		interp_out = @(x) NaN;
	end
end

function interp_out = constrained_interp(values, cdf_x)
	[values_u, iu] = unique(values, 'last');
	cdf_x_u = cdf_x(iu);

	if numel(values_u) >= 2
		interp_out = griddedInterpolant(...
			values_u, cdf_x_u, 'pchip', 'nearest');
	else
		interp_out = @(x) NaN;
	end
end

% function [values_out, cdf_out, iu] = unique_sort(values_in, pmf_in, iunique)
% 	sorted_mat = sortrows([values_in(:), pmf_in(:)]);
% 	tmp_cdf = cumsum(sorted_mat(:,2));

% 	if iunique == 1
% 		[values_out, iu] = unique(sorted_mat, 'last');
% 		cdf_out = tmp_cdf(iu);
% 	else
% 		[cdf_out, iu] = unique(tmp_cdf, 'first');
% 		values_out = sorted_mat(iu,1);
% 	end
% end