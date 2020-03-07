classdef TableGen < handle

	properties (SetAccess = protected)
		mpcs_present = false;
		mpcs_news_present = false;
		mpcs_loan_loss_present = false;
		decomp_norisk_present = false;
		decomp_RA_present = false;
		selected_cases;
		freq;
	end

	properties
		outdir;
		output = table();
		decomp_baseline;
		decomp_incrisk_alt;
	end

	properties (Abstract)
		default_fname;
		included_names;
	end

	methods
		function obj = TableGen(params, results, freq, use_all)
			obj.freq = freq;

			if ~isequal(freq, [1, 4]);
				this_freq = find([params.freq] == freq);
				if nargin < 4
					use_all = true;
				end
			else
				this_freq = 1:numel(params);
			end

			if use_all || isempty(obj.included_names)
				obj.selected_cases = this_freq;
			else
				all_names = [params.name];
				inames = [];
				for ii = 1:numel(obj.included_names)
					name = obj.included_names(ii);
					inames = [inames, find(ismember(all_names, name))];
				end

				obj.selected_cases = intersect(this_freq, inames);
			end

			obj.check_options(params, results);
			obj.outdir = params(1).outdir;
		end

		function check_options(obj, params, results)

			if isempty(obj.selected_cases)
				return
			end

			for ip = obj.selected_cases
				if params(ip).MPCs
					obj.mpcs_present = true;
				end

				if params(ip).MPCs_news
					obj.mpcs_news_present = true;
				end

				if params(ip).MPCs_loan_and_loss
					obj.mpcs_loan_loss_present = true;
				end

				if results(ip).decomp_norisk(1).completed
					obj.decomp_norisk_present = true;
				end

				if results(ip).decomp_RA.completed
					obj.decomp_RA_present = true;
				end
			end
		end

		function save_table(obj, fname)
			if nargin == 1
				fpath = fullfile(obj.outdir, obj.default_fname);
			else
				fpath = fullfile(obj.outdir, fname);
			end
			writetable(obj.output, fpath, 'WriteRowNames', true,...
				'WriteVariableNames', false);
		end
	end

	methods (Static)
		function output_table = append_to_table(...
				input_table, new_entries, new_labels)
			table_to_append = table(new_entries,...
				'VariableNames', {'results'},...
				'RowNames', new_labels);
			output_table = [input_table; table_to_append];
		end

		function new_table = new_table_with_header(header_name)
			header_formatted = strcat('____', header_name);
			new_table = table({NaN},...
				'VariableNames', {'results'},...
				'RowNames', {header_formatted});
		end
	end
end