classdef Calibrator
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	methods (Static)
		function x = fraction_constrained(beta_in, p)
			
			p.set("beta0", beta_in);
			results = main(p);

			x = results.direct.wealth_lt_1000 - 0.23;

			fprintf("\n\n---- For beta = %f, P(assets<$1000) = %f ----\n\n",...
				p.beta0, results.direct.wealth_lt_1000)

		end

		function x = mean_wealth(input_value, p, input_name)
			p.set(input_name, input_value);

			if strcmp(input_name, "r")
				p.set("R", 1 + input_value);
			end

			results = main(p);

			x = results.direct.mean_a - 0.25;

			fprintf("\n\n---- For %s = %f, E[a] = %f ----\n\n",...
				input_name, input_value, results.direct.mean_a)
		end
	end
end
