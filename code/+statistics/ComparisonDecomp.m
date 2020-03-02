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

		agrid;
		incompatible_params = false;

		Empc0_g0;
		Empc0_g1;
		Empc1_g0;
		Empc1_g1;

		integral_m0g0_interp;
		integral_m0g1_interp;
		integral_m1g0_interp;
		integral_m1g1_interp;

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

			obj.agrid = stats0.agrid;
			obj.na = numel(obj.agrid);

			obj.nthresholds = numel(obj.p0.abars);

			obj.initialize();
		end

		function initialize(obj)
			obj.results.Em1_less_Em0 = NaN;
			obj.results.term1 = NaN;
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

			obj.Empc0_g0 = dot(mpcs0_a, obj.pmf0_a);
			obj.Empc0_g1 = dot(mpcs0_a, obj.pmf1_a);
			obj.Empc1_g0 = dot(mpcs1_a, obj.pmf0_a);
            obj.Empc1_g1 = dot(mpcs1_a, obj.pmf1_a);

            obj.integral_m0g0_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs0_a, obj.pmf0_a, true);
            obj.integral_m0g1_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs0_a, obj.pmf1_a, true);
            obj.integral_m1g0_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs1_a, obj.pmf0_a, true);
            obj.integral_m1g1_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs1_a, obj.pmf1_a, true);
		end
	end

end