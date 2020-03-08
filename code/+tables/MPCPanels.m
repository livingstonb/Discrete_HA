classdef MPCPanels
	methods (Static)
		function out = size_effects(values, shocks_labels, include_annual,...
			panel_name)
			if nargin < 4
				panel_name = 'Panel C: MPC Size Effects';
			end

			out = tables.TableGen.new_table_with_header(panel_name);

			if include_annual
				new_labels = {	sprintf('Annual MPC (%%), %s', shocks_labels{4})
				            	sprintf('Annual MPC (%%), %s', shocks_labels{6})
					};
				new_entries = {	values.direct.mpcs(4).avg_annual * 100
								values.direct.mpcs(6).avg_annual * 100
					};
				new_entries = aux.cellround(new_entries, 1);
				out = tables.TableGen.append_to_table(out, new_entries, new_labels);
			end

			new_labels = {	sprintf('Quarterly MPC (%%), %s', shocks_labels{4})
				            sprintf('Quarterly MPC (%%), %s', shocks_labels{6})
				};
			new_entries = {	values.direct.mpcs(4).avg_quarterly * 100
							values.direct.mpcs(6).avg_quarterly * 100
				};
			new_entries = aux.cellround(new_entries, 1);
			out = tables.TableGen.append_to_table(out, new_entries, new_labels);
		end

		function out = sign_effects(values, shocks_labels, include_annual,...
			panel_name)
			if nargin < 4
				panel_name = 'Panel D: MPC Sign Effects';
			end

			out = tables.TableGen.new_table_with_header(panel_name);

			if include_annual
				new_labels = {	sprintf('Annual MPC (%%), %s', shocks_labels{1})
								sprintf('Annual MPC (%%), %s', shocks_labels{2})
					            sprintf('Annual MPC (%%), %s', shocks_labels{3})
					};
				new_entries = {	values.direct.mpcs(1).avg_annual * 100
								values.direct.mpcs(2).avg_annual * 100
								values.direct.mpcs(3).avg_annual * 100
					};
				new_entries = aux.cellround(new_entries, 1);
				out = tables.TableGen.append_to_table(out, new_entries, new_labels);
			end

			new_labels = {	sprintf('Quarterly MPC (%%), %s', shocks_labels{1})
				            sprintf('Quarterly MPC (%%), %s', shocks_labels{2})
				            sprintf('Quarterly MPC (%%), %s', shocks_labels{3})
				};
			new_entries = {	values.direct.mpcs(1).avg_quarterly * 100
							values.direct.mpcs(2).avg_quarterly * 100
							values.direct.mpcs(3).avg_quarterly * 100
				};
			new_entries = aux.cellround(new_entries, 1);
			out = tables.TableGen.append_to_table(out, new_entries, new_labels);
		end
	end
end