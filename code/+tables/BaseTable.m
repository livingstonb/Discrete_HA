classdef BaseTable < handle
	properties (SetAccess = protected)
		mpcs_present = false;
		mpcs_news_present = false;
		mpcs_loan_loss_present = false;
		decomp_norisk_present = false;
		decomp_RA_present = false;
		
		n_cols;
		current_column;
	end

	properties
		outdir;
		output;
		decomp_incrisk_alt;
		decomp_baseline;
	end

	methods
		function obj = BaseTable(params, stats)
            obj.outdir = params(1).outdir;
			obj.n_cols = numel(params);

			obj.set_options(params, stats);
		end

		function set_options(obj, params, stats)
			obj.mpcs_present = any([params.MPCs]);
			obj.mpcs_news_present = any([params.MPCs_news]);
			obj.mpcs_loan_loss_present = any([params.MPCs_loan_and_loss]);

% 			for ip = 1:numel(params)
% 				if stats{ip}.decomp_norisk(1).completed
% 					obj.decomp_norisk_present = true;
% 				end
% 
% 				if stats{ip}.decomp_RA.completed
% 					obj.decomp_RA_present = true;
% 				end
% 			end
		end

		function update_current_column(obj, table_in, stats_in)
			if nargin < 3
				tmp = table_in;
			else
				tmp = obj.construct_from_stats(table_in, stats_in);
			end
			obj.current_column = [obj.current_column; tmp];
		end

		function table_out = construct_from_stats(obj, table_in, stats_in)
			vals = {};
			labels = {};
			for ii = 1:numel(stats_in)
				vals{ii} = stats_in{ii}.value;
				labels{ii} = stats_in{ii}.label;
			end

			table_to_append = table(vals(:),...
				'VariableNames', {'results'},...
				'RowNames', labels(:));

			table_out = [table_in; table_to_append];
		end

		function add_column(obj, ip)
			column_label = sprintf('Specification%d', ip);
			obj.current_column.Properties.VariableNames = {column_label};
			obj.output = [obj.output, obj.current_column];
		end

		function new_table = new_table_with_header(obj, header_name)
			header_formatted = strcat('____', header_name);
			new_table = table({NaN},...
				'VariableNames', {'results'},...
				'RowNames', {header_formatted});
		end
	end

	methods (Static)
		function s_out = sround(sstruct, n)
			s_out = sstruct;
			s_out.value = round(s_out.value, n);
		end

		function s_out = sround_mult(scell, n)
			s_out = cell(numel(scell), 1);
			for ii = 1:numel(scell)
				s_out{ii} = scell{ii};
				s_out{ii}.value = round(s_out{ii}.value, n);
			end
		end
	end
end