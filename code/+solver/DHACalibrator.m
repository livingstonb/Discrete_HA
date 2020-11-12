classdef DHACalibrator < EconToolsML.Calibrator
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	methods
		function obj = DHACalibrator(params, variables,...
			target_names, target_values)

			obj = obj@EconToolsML.Calibrator(params, variables, target_names, target_values);
			
			obj.main_handle = @(curr_params) main(curr_params, 'iterating', true);
		end

		function construct_options_struct(obj, params)
			obj.options = struct();
			obj.options.MPCs = params.MPCs;
			obj.options.MPCs_news = params.MPCs_news;
			obj.options.MPCs_loan_and_loss = params.MPCs_loan_and_loss;
			obj.options.Simulate = params.Simulate;
			obj.options.DeterministicMPCs = params.DeterministicMPCs;
			obj.options.MakePlots = params.MakePlots;
			obj.options.SaveOutput = params.SaveOutput;
		end

		function value = get_results_value(obj, results, variable_name)
			value = results.direct.(variable_name);
		end
	end
end