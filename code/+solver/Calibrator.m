classdef Calibrator < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties
		options;
		variable;
		target_name;
		target_value;
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
			obj.variable = variable;

			obj.target_name = target_name;
			obj.target_value = target_value;
		end

		function dv = fn_handle(obj, x, current_params)
			obj.turn_off_param_options(current_params);
			current_params.set(obj.variable, x);
			results = main(current_params);

			v = results.direct.(obj.target_name);
			dv = v - obj.target_value;

			fprintf('\n\n---- For %s = %g:\n', obj.variable, x)
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
	end
end
