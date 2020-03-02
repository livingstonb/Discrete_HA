classdef Decomp < handle
	properties
		agrid;
		stats;

		na;
		nthresholds;

		pmf;
		pmf_a;

		integral_interp;
		integral_norisk_interp;
		cdf_interp;

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
		end

		function make_initial_computations(obj, mpcs)
			obj.results_RA = struct('mpc_less_mpcRA', NaN,...
				'term1', NaN, 'term2', NaN, 'term3');

			tmp = (1-obj.p.dieprob) * obj.stats.beta * obj.p.R;
            mpc_ra = obj.p.R * tmp ^ (-1 / obj.p.risk_aver) - 1;

            obj.pmf = obj.stats.adist;
            obj.pmf_a = obj.stats.agrid_dist;

            mpcs_a = obj.collapse_mpcs(mpcs)

            obj.integral_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs_a, obj.pmf_a, true);

            obj.integral_norisk_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs_norisk, obj.pmf_a, true);

            cdf_a = cumsum(obj.pmf_a);
            dsupport = obj.pmf_a > 1e-7;
            obj.cdf_interp = griddedInterpolant(...
            	obj.agrid(dsupport), cdf_a(dsupport), 'linear');
		end

		function perform_decompositions(obj, mpcs, mpcs_norisk)
			obj.make_initial_computations(mpcs, mpcs_norisk);
			obj.decomp_RA();
			obj.decomp_norisk();
		end

		function decomp_RA(obj)
		end

		function decomp_norisk(obj)
		end

		function mpcs_a = collapse_mpcs(obj, mpcs_states)
			mpcs_states = reshape(mpcs_states, obj.na, []);
			mpcs_a = sum(mpcs_states .* obj.pmf, 2)...
				./ obj.pmf_a;

			pmf_a_small = obj.pmf_a < 1e-8;
			mpcs_a(pmf_a_small) = mean(mpcs_states(pmf_a_small,:), 2);
		end
	end
end