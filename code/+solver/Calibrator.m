classdef Calibrator < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties
		options;
		variables;
		target_names;
		target_values;
		lbounds = [];
		ubounds = [];
		n;
		x0;

		solver_handle;
	end

	methods
		function obj = Calibrator(params, variables,...
			target_names, target_values);
			options = struct();
			obj.options.MPCs = params.MPCs;
			obj.options.MPCs_news = params.MPCs_news;
			obj.options.MPCs_loan_and_loss = params.MPCs_loan_and_loss;
			obj.options.Simulate = params.Simulate;
			obj.options.DeterministicMPCs = params.DeterministicMPCs;
			obj.options.MakePlots = params.MakePlots;

			obj.variables = variables;
			obj.n = numel(variables);

			obj.target_names = target_names;
			obj.target_values = target_values;

			if obj.n ~= numel(obj.target_names)
				error("Number of instruments and targets don't match")
			elseif numel(obj.target_values) ~= numel(obj.target_names)
				error("Too many/few values provided for targets")
			end

			for i_var = 1:obj.n
				obj.x0(i_var) = params.(obj.variables{i_var});
			end
				
		end

		function set_param_bounds(obj, varargin)
			nv = numel(varargin);
			assert(nv==obj.n, "Too many or too few bounds provided");

			for ii = 1:nv
				var_bounds = varargin{ii};
				obj.lbounds(ii) = var_bounds(1);
				obj.ubounds(ii) = var_bounds(2);
			end
		end

		function set_handle(obj, p)
			for ii = 1:obj.n
				if obj.x0(ii) < obj.lbounds(ii)
					obj.x0(ii) = (5*obj.lbounds(ii) + obj.ubounds(ii))/6;
				elseif obj.x0(ii) > obj.ubounds(ii)
					obj.x0(ii) = (obj.lbounds(ii) + 5*obj.ubounds(ii))/6;
				end
			end

			obj.solver_handle = @(x) obj.fn_handle(x, p);
		end

		function dv = fn_handle(obj, x, current_params)
			obj.turn_off_param_options(current_params);

			for i_var = 1:obj.n
				current_params.set(obj.variables{i_var}, x(i_var));
			end
			results = main(current_params);

			fprintf('\n\n---- For ')
			for i_var = 1:obj.n
				v(i_var) = results.direct.(obj.target_names{i_var});
				fprintf('%s = %g', obj.variables{i_var}, x(i_var))
				if i_var < obj.n
					fprintf(", ")
				else
					fprintf(':\n')
				end
			end
			dv = v - obj.target_values;

			fprintf('    Result was ')
			for i_var = 1:obj.n
				fprintf('%s = %g', obj.target_names{i_var}, v(i_var))
				if i_var < obj.n
					fprintf(", ")
				else
					fprintf('\n')
				end
			end

			obj.reset_param_options(current_params);
		end

		function solver_args = get_args(obj)
			if ~isempty(obj.ubounds)
				solver_args = {obj.x0, obj.lbounds, obj.ubounds};
			elseif ~isempty(obj.x0)
				solver_args = {obj.x0};
			else
				solver_args = {};
			end
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
