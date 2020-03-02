classdef ComparisonDecomp

	properties
		p0;
		p1;
		stats0;
		stats1;

		pmf0;
		pmf1;
		pmf0_a;
		pmf1_a;

		can_compute_RAmpc;
		mpc0_ra = NaN;
		mpc1_ra = NaN;

		agrid;
		incompatible_params = false;

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

			if ~isequal(stats0.agrid, stats1.agrid)
				obj.incompatible_params = true;
			elseif ~isequal(p0.abars, p1.abars)
				obj.incompatible_params = true;
			end

			no_z_heterogeneity = all([p0.nb, p1.nb] == 1);
			no_temptation = isequal(p0.temptation, 0)...
				&& isequal(p1.temptation, 0);
			no_bequest_prefs = isequal(p0.bequest_weight, 0)...
				&& isequal(p1.bequest_weight, 0);
			no_ez = all(~[p0.EpsteinZin, p1.EpsteinZin]);

			obj.can_compute_RAmpc = no_z_heterogeneity...
				&& no_temptation && no_bequest_prefs...
				&& no_ez;

			obj.agrid = stats0.agrid;
			obj.na = numel(obj.agrid);

			obj.nthresholds = numel(obj.p0.abars);

			obj.initialize();
		end

		function initialize(obj)
			obj.results.Em1_less_Em0 = NaN;
			obj.results.term1 = NaN;
			obj.results.term1a = NaN;
			obj.results.term1b = NaN;
			obj.results.term2 = NaN;
			obj.results.term2a = NaN(1, obj.nthresholds);
			obj.results.term2b = NaN(1, obj.nthresholds);
			obj.results.term3 = NaN;
			obj.results.complete = false;
		end

		function perform_decompositions(obj, mpcs0, mpcs1)
			if obj.incompatible_params
				return
			end

			obj.make_initial_computations(mpcs0, mpcs1);

			obj.results.Em1_less_Em0 = obj.Empc1_g1 - obj.Empc1_g0;

			% Term 1: Effect of MPC function
			obj.results.term1 = obj.Empc1_g1 - obj.Empc0_g1;

			if obj.can_compute_RAmpc
				offset = obj.mpc1_ra - obj.mpc0_ra;


				% Term 1a: Effect of MPC function, level
				obj.results.term1a = obj.Empc1_g0 - obj.Empc1adj_g0;

				% Term 1b: Effect of MPC function, shape
				obj.results.term1b = obj.Empc1adj_g0 - obj.Empc0_g0;
			end


			% Term 2: Effect of distribution
			obj.results.term2 = obj.Empc0_g1 - obj.Empc0_g0;

			% Term 3: Interaction
			obj.results.term2 = (obj.Empc1_g1 - obj.Empc0_g1)...
				- (obj.Empc1_g0 - obj.Empc0_g0);


			for ia = 1:obj.nthresholds
				thresh = obj.p.abars(ia);

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
			obj.pmf = obj.stats.adist;
            obj.pmf_a = obj.stats.agrid_dist;

			obj.pmf0 = obj.stats0.adist;
			obj.pmf0_a = obj.stats0.agrid_dist;
			obj.pmf1 = obj.stats1.adist;
			obj.pmf1_a = obj.stats1.agrid_dist;

			mpcs0_a = aux.collapse_mpcs(mpcs0, obj.pmf0, obj.pmf0_a);
			mpcs1_a = aux.collapse_mpcs(mpcs1, obj.pmf1, obj.pmf1_a);

			if obj.can_compute_RAmpc
				tmp0 = (1-obj.p0.dieprob) * obj.stats0.beta * obj.p0.R;
	            obj.mpc0_ra = obj.p0.R * tmp ^ (-1 / obj.p0.risk_aver) - 1;

	            tmp1 = (1-obj.p1.dieprob) * obj.stats1.beta * obj.p1.R;
	            obj.mpc1_ra = obj.p1.R * tmp ^ (-1 / obj.p1.risk_aver) - 1;
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