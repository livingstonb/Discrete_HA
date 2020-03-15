classdef TableGen < handle

	properties (SetAccess = protected)
		mpcs_present = false;
		mpcs_news_present = false;
		mpcs_loan_loss_present = false;
		decomp_norisk_present = false;
		decomp_RA_present = false;
		this_freq;
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
	end

	methods
		function check_options(obj, params, results)

			for ip = obj.this_freq
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