classdef ComparisonDecomp < handle

	properties
		na;
		nthresholds;
		thresholds;

		p0;
		p1;
		stats0;
		stats1;

		pmf0_a;
		pmf1_a;

		can_compute_RAmpc;
		mpc0_ra = NaN;
		mpc1_ra = NaN;

		agrid;
		incompatible_params = false;
		same_model_name;

		Empc0_g0;
		Empc0_g1;
		Empc1_g0;
		Empc1_g1;
		Empc1adj_g0 = NaN;

		integral_m0g0_interp;
		integral_m0g1_interp;

		results = struct();
	end

	methods
		function obj = ComparisonDecomp(p0, p1, stats0, stats1)
			obj.p0 = p0;
			obj.p1 = p1;
			obj.stats0 = stats0;
			obj.stats1 = stats1;

			if ~isequal(p0.nx_DST, p1.nx_DST)
				obj.incompatible_params = true;
			elseif ~isequal(p0.abars, p1.abars)
				obj.incompatible_params = true;
			end

			obj.same_model_name = strcmp(p0.name, p1.name);

			no_z_heterogeneity = all([p0.nb, p1.nb] == 1);
			no_temptation = isequal(p0.temptation, 0)...
				&& isequal(p1.temptation, 0);
			no_bequest_prefs = isequal(p0.bequest_weight, 0)...
				&& isequal(p1.bequest_weight, 0);
			no_ez = all([~p0.EpsteinZin, ~p1.EpsteinZin]);

			obj.can_compute_RAmpc = no_z_heterogeneity...
				&& no_temptation && no_bequest_prefs...
				&& no_ez;

			obj.agrid = stats0.agrid;
			obj.na = numel(obj.agrid);

			obj.nthresholds = numel(obj.p0.abars);
			obj.thresholds = obj.p0.abars;

			obj.initialize();
		end

		function initialize(obj)
			obj.results.Em1_less_Em0 = NaN;
			obj.results.term1 = NaN;
			obj.results.term1a = NaN;
			obj.results.term1b = NaN;
			obj.results.term2 = NaN;
			obj.results.term2a = NaN(obj.nthresholds, 1);
			obj.results.term2b = NaN(obj.nthresholds, 1);
			obj.results.term3 = NaN;
			obj.results.complete = false;
		end

		function perform_decompositions(obj, mpcs0, mpcs1)
			if obj.incompatible_params || obj.same_model_name
				return
			end

			if ~all([obj.p0.MPCs, obj.p1.MPCs,...
					obj.p0.DeterministicMPCs, obj.p1.DeterministicMPCs]);
				return
			end

			obj.make_initial_computations(mpcs0, mpcs1);

			obj.results.Em1_less_Em0 = obj.Empc1_g1 - obj.Empc0_g0;

			% Term 1: Effect of MPC function
			obj.results.term1 = obj.Empc1_g0 - obj.Empc0_g0;

			if obj.can_compute_RAmpc
				% Term 1a: Effect of MPC function, level
				obj.results.term1a = obj.Empc1_g0 - obj.Empc1adj_g0;

				% Term 1b: Effect of MPC function, shape
				obj.results.term1b = obj.Empc1adj_g0 - obj.Empc0_g0;
			end

			% Term 2: Effect of distribution
			obj.results.term2 = obj.Empc0_g1 - obj.Empc0_g0;

			% Term 3: Interaction
			obj.results.term3 = (obj.Empc1_g1 - obj.Empc0_g1)...
				- (obj.Empc1_g0 - obj.Empc0_g0);

			for ia = 1:obj.nthresholds
				thresh = obj.thresholds(ia);

				% Term 2a: Dist effect for HtM households
				obj.results.term2a(ia) = obj.integral_m0g1_interp(thresh)...
					- obj.integral_m0g0_interp(thresh);

				% Term 2b: Dist effect for NHtM households
				obj.results.term2b(ia) =...
					(obj.Empc0_g1 - obj.integral_m0g1_interp(thresh))...
					- (obj.Empc0_g0 - obj.integral_m0g0_interp(thresh));
			end

			obj.results.complete = true;
		end

		function make_initial_computations(obj, mpcs0, mpcs1)

			pmf0 = obj.stats0.adist;
			obj.pmf0_a = obj.stats0.agrid_dist;
			mpcs0_a = aux.collapse_mpcs(mpcs0, pmf0, obj.pmf0_a);

			if ~isequal(obj.stats0.agrid, obj.stats1.agrid)
				% Need to reconstruct mpc's and pmf for model 1 onto
				% baseline grid

				agrid1_orig = obj.stats1.agrid;
				pmf1_orig = obj.stats1.adist;
				pmf1_a_orig = obj.stats1.agrid_dist;
				mpcs1_a_orig = aux.collapse_mpcs(...
					mpcs1, pmf1_orig, pmf1_a_orig);

				% First get cdf interpolant
                cdf1_a_interp = griddedInterpolant(...
                	agrid1_orig, cumsum(pmf1_a_orig), 'pchip', 'nearest');

                cdf_a0 = pmf1_a_orig(1);
                cdf1_a_interp = @(x) adjust_interpolant(x,...
                	cdf1_a_interp, agrid1_orig, cdf_a0);

            	% Next get pmf on the baseline grid
            	obj.pmf1_a = zeros(obj.na, 1);
            	obj.pmf1_a(1) = cdf1_a_interp(obj.agrid(1));
            	for ia = 2:obj.na
            		obj.pmf1_a(ia) = cdf1_a_interp(obj.agrid(ia)) - cdf1_a_interp(obj.agrid(ia-1));
            	end
            	obj.pmf1_a = obj.pmf1_a / sum(obj.pmf1_a);

            	% Next get mpc function on baseline grid
            	mpcs1_a_interp = griddedInterpolant(...
            		agrid1_orig, mpcs1_a_orig, 'pchip', 'nearest');
            	mpcs1_a = mpcs1_a_interp(obj.agrid);
            else
            	pmf1 = obj.stats1.adist;
				obj.pmf1_a = obj.stats1.agrid_dist;
				mpcs1_a = aux.collapse_mpcs(mpcs1, pmf1, obj.pmf1_a);
			end

			if obj.can_compute_RAmpc
				tmp0 = (1-obj.p0.dieprob) * obj.stats0.beta * obj.p0.R;
	            obj.mpc0_ra = obj.p0.R * tmp0 ^ (-1 / obj.p0.risk_aver) - 1;

	            tmp1 = (1-obj.p1.dieprob) * obj.stats1.beta * obj.p1.R;
	            obj.mpc1_ra = obj.p1.R * tmp1 ^ (-1 / obj.p1.risk_aver) - 1;
	            offset = obj.mpc1_ra - obj.mpc0_ra;

	            mpcs1_adj = mpcs1_a - offset;

	            obj.Empc1adj_g0 = dot(mpcs1_adj, obj.pmf0_a);
	        end

			obj.Empc0_g0 = dot(mpcs0_a, obj.pmf0_a);
			obj.Empc0_g1 = dot(mpcs0_a, obj.pmf1_a);
			obj.Empc1_g0 = dot(mpcs1_a, obj.pmf0_a);
            obj.Empc1_g1 = dot(mpcs1_a, obj.pmf1_a);

            obj.integral_m0g0_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs0_a, obj.pmf0_a, true);
            obj.integral_m0g1_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs0_a, obj.pmf1_a, true);
		end
	end
end

function vals_out = adjust_interpolant(x, cdf1_a_interp, agrid1_orig, cdf_a0)
	vals_out = cdf1_a_interp(x);

	x0 = [0; agrid1_orig(1)];
	vals0 = [0; cdf_a0];

	low_states = x < agrid1_orig(1);
    
    if sum(low_states) > 0
        vals_out(low_states) = interp1(x0, vals0, x(low_states));
    end
end