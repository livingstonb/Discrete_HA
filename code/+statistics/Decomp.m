classdef Decomp < handle
	properties
		agrid;
		stats;
        p;

		na;
		nthresholds;

		pmf;
		pmf_a;

		integral_interp;
		integral_norisk_interp;
		cdf_interp;

		Empc;
		Empc_norisk;
		mpc_ra;
		mpcs_a;
		mpcs_norisk;

		results_RA;
		results_norisk;
	end

	methods
		function obj = Decomp(p, stats)
			obj.p = p;
			obj.agrid = stats.agrid;
			obj.na = numel(obj.agrid);
			obj.stats = stats;

			obj.nthresholds = numel(obj.p.abars);

			obj.initialize();
		end

		function perform_decompositions(obj, mpcs, mpcs_norisk)
			obj.mpcs_norisk = mpcs_norisk;

			obj.make_initial_computations(mpcs);
			obj.decomp_RA();
			obj.decomp_norisk();
		end

		function initialize(obj)
			obj.results_norisk = struct();
			for ia = 1:obj.nthresholds
				obj.results_norisk(ia).completed = false;
				obj.results_norisk(ia).term1 = NaN;
				obj.results_norisk(ia).term2 = NaN;
				obj.results_norisk(ia).term3 = NaN;
				obj.results_norisk(ia).term4 = NaN;
			end

			obj.results_RA = struct();
			obj.results_RA.completed = false;
			obj.results_RA.Em1_less_mRA = NaN;
			obj.results_RA.term1 = NaN;
			obj.results_RA.term2 = NaN;
			obj.results_RA.term3 = NaN;
		end

		function make_initial_computations(obj, mpcs)
			tmp = (1-obj.p.dieprob) * obj.stats.beta * obj.p.R;
            obj.mpc_ra = obj.p.R * tmp ^ (-1 / obj.p.risk_aver) - 1;

            obj.pmf = obj.stats.adist;
            obj.pmf_a = obj.stats.agrid_dist;

            obj.mpcs_a = aux.collapse_mpcs(mpcs, obj.pmf, obj.pmf_a);

            obj.Empc = dot(obj.mpcs_a, obj.pmf_a);
            obj.Empc_norisk = dot(obj.mpcs_norisk, obj.pmf_a);

            obj.integral_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs_a, obj.pmf_a, true);

            obj.integral_norisk_interp = aux.interpolate_integral(...
            	obj.agrid, obj.mpcs_norisk, obj.pmf_a, true);

            integrand = ones(size(obj.pmf_a));
            obj.cdf_interp = aux.interpolate_integral(...
            	obj.agrid, integrand, obj.pmf_a, true);
		end

		function decomp_RA(obj)
			obj.results_RA.Em1_less_mRA = obj.Empc - obj.mpc_ra;

			% Interpolant to compute E[MPC|a=E[a]]
			dsupport = obj.pmf_a > 1e-7;
			mpc_a_interp = griddedInterpolant(...
				obj.agrid(dsupport), obj.mpcs_a(dsupport), 'linear');
			mpc_atmean = mpc_a_interp(obj.stats.mean_a);

			% Term 1: Effect of MPC function
			obj.results_RA.term1 = mpc_atmean - obj.mpc_ra;

			% Term 2: Effect of distribution
			obj.results_RA.term2 = 0;

			% Term 3: Interaction
			obj.results_RA.term3 = (obj.Empc - obj.mpc_ra)...
				- (mpc_atmean - obj.mpc_ra);

			obj.results_RA.completed = true;
		end

		function decomp_norisk(obj)
			for ia = 1:obj.nthresholds
				% Term 1: RA mpc
				obj.results_norisk(ia).term1 = obj.mpc_ra;

				thresh = obj.p.abars(ia);
				if thresh == min(obj.agrid)
					% Term 2: HtM effect
					obj.results_norisk(ia).term2 = ...
						(obj.mpcs_a(1) - obj.mpc_ra) * obj.pmf_a(1);

					% Term 3: Constraint effect for non-HtM
					obj.results_norisk(ia).term3 = ...
						dot(obj.mpcs_norisk(2:end) - obj.mpc_ra, obj.pmf_a(2:end));

					% Term 4: Income risk effect for non-HtM
					obj.results_norisk(ia).term4 = ...
						dot(obj.mpcs_a(2:end), obj.pmf_a(2:end))...
						- dot(obj.mpcs_norisk(2:end), obj.pmf_a(2:end));
				else
					% Term 2: HtM effect
					obj.results_norisk(ia).term2 = ...
						obj.integral_interp(thresh) - obj.mpc_ra * obj.cdf_interp(thresh);

					% Term 3: Constraint effect for non-HtM
					obj.results_norisk(ia).term3 = ...
						(obj.Empc_norisk - obj.integral_norisk_interp(thresh))...
						- obj.mpc_ra * (1 - obj.cdf_interp(thresh));

					% Term 4: Income risk effect for non-HtM
					obj.results_norisk(ia).term4 = ...
						(obj.Empc - obj.integral_interp(thresh))...
						- (obj.Empc_norisk - obj.integral_norisk_interp(thresh));
				end
				obj.results_norisk(ia).completed = true;
			end
		end
	end
end