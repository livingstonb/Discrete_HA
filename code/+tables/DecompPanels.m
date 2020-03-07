classdef DecompPanels
	methods (Static)
		function out = decomp_intro(values, params)
			out = table({params.name},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_labels = {	'MPC, quarterly or annual (%)'
							'RA MPC (%)'
				};
			new_entries = {	values.mpcs(5).avg_s_t(1,1) * 100
							values.mpc_RA * 100
				};

			new_entries = aux.cellround(new_entries, 1);
			out = tables.TableGen.append_to_table(...
				out, new_entries, new_labels);
		end

		function out = mean_mpc_decomp(decomp, params, ithresh, panel_prefix)
			thresh = params.convert_to_dollars_str(params.abars(ithresh));

			panel_suffix = sprintf(...
				'Decomposition of Mean MPC, around %s', thresh);
			if nargin < 4
				panel_name = panel_suffix;
			else
				panel_name = strcat(panel_prefix, ": ", panel_suffix);
			end
			
			out = tables.TableGen.new_table_with_header(char(panel_name));
			new_labels = {	sprintf('HtM (a <= %s), effect ', thresh)
				            sprintf('Non-HtM (a > %s), constraint effect', thresh)
				            sprintf('Non-HtM (a > %s), inc risk effect', thresh)
				};
			new_entries = {	decomp(ithresh).term2
							decomp(ithresh).term3
							decomp(ithresh).term4
				};

			new_entries = aux.cellround(new_entries, 3);
			out = tables.TableGen.append_to_table(...
				out, new_entries, new_labels);
		end

		function out = ra_mpc_decomp(decomp, panel_prefix)
			panel_suffix = 'Decomposition of Mean MPC - MPC_RA';
			if nargin < 2
				panel_name = panel_suffix;
			else
				panel_name = strcat(panel_prefix, ": ", panel_suffix);
			end

			out = tables.TableGen.new_table_with_header(char(panel_name));
			new_labels = {	'Mean MPC - RA MPC'
				            'Effect of MPC function'
				            'Effect of distribution'
				            'Interaction'
				};
			new_entries = {	decomp.Em1_less_mRA
							decomp.term1
		                    decomp.term2
		                    decomp.term3
				};
			new_entries = aux.cellround(new_entries, 3);
			out = tables.TableGen.append_to_table(...
				out, new_entries, new_labels);
		end
	end
end