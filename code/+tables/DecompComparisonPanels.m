classdef DecompComparisonPanels
	methods (Static)
		function out = decomp_wrt_baseline(values, params, ishock,...
			absolute, panel_prefix)

			shock_size = params.shocks(ishock);
			shock_str = params.convert_to_dollars_str(shock_size);

			if absolute
				panel_suffix = 'Quarterly MPC decomposition relative to baseline';
			else
				panel_suffix = ...
					'Quarterly MPC decomposition, % of E[MPC] - E[MPC_baseline]';
			end

			if nargin < 5
				panel_name = panel_suffix;
			else
				panel_name = strcat(panel_prefix, ": ", panel_suffix);
			end

			out = tables.TableGen.new_table_with_header(char(panel_name));

			if absolute
				new_labels = {'E[MPC] - E[MPC_baseline]'};
				new_entries = {values.Em1_less_Em0};
				new_entries = aux.cellround(new_entries, 3);
				out = tables.TableGen.append_to_table(...
					out, new_entries, new_labels);
			end

			new_labels = {	'Effect of MPC fcn'
							'Effect of distr'
							'  Distr effect, HtM'
							'  Distr effect, NHtM'
							'Interaction'
				};
			new_entries = [	values.term1
							values.term2
							values.term2a(2)
							values.term2b(2)
							values.term3
				];

			if absolute
				new_entries = aux.cellround(...
					num2cell(new_entries), 3);
			else
				new_labels = cellfun( @(entry) ...
					char(strcat(entry, " (%)")), new_labels, 'UniformOutput', false);

				new_entries = new_entries / values.Em1_less_Em0 * 100;
				new_entries = aux.cellround(num2cell(new_entries), 1);
			end

			out = tables.TableGen.append_to_table(...
				out, new_entries, new_labels);
		end
	end
end
