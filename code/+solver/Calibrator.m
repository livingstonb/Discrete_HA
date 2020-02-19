classdef Calibrator < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties
		options;
		variable;
		target_name;
		target_value;
		lbound;
		ubound;
	end

	methods
		function obj = Calibrator(params, variable,...
			target_name, target_value);
			options = struct();
			obj.options.MPCs = params.MPCs;
			obj.options.MPCs_news = params.MPCs_news;
			obj.options.MPCs_loan_and_loss = params.MPCs_loan_and_loss;
			obj.options.Simulate = params.Simulate;
			obj.options.DeterministicMPCs = params.DeterministicMPCs;
			obj.options.MakePlots = params.MakePlots;

			obj.variable = variable;

			obj.target_name = target_name;
			obj.target_value = target_value;
		end

		function set_param_bounds(obj, bounds)
			obj.lbound = bounds(1);
			obj.ubound = bounds(2);
		end

		function dv = fn_handle(obj, x, current_params)
			obj.turn_off_param_options(current_params);

			if ~isempty(obj.lbound)
				new_value = obj.convert_to_param_value(x);
			else
				new_value = x;
			end

			current_params.set(obj.variable, new_value);
			results = main(current_params);

			v = results.direct.(obj.target_name);
			dv = v - obj.target_value;

			fprintf('\n\n---- For %s = %g:\n', obj.variable, new_value)
			fprintf('     Result was %s = %g\n', obj.target_name, v)

			obj.reset_param_options(current_params);
		end

		function turn_off_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, 0, quiet);
			end
		end

		function reset_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, obj.options.(props{ip}), quiet);
			end
		end

		function val = convert_to_param_value(obj, z)
			lb = obj.lbound;
			ub = obj.ubound;
			tmp = abs(z) / (1 + abs(z));
			val = lb + tmp * (ub - lb);
		end

		function z = convert_to_solver_input(obj, param_value)
			lb = obj.lbound;
			ub = obj.ubound;
			tmp = (param_value - lb) / (ub - lb);
			z = tmp / (1 - tmp);
		end
	end
end
