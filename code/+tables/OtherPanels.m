classdef OtherPanels
	methods (Static)
		function out = intro_panel(values, p)
			out = table({p.name},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_labels = {	'Quarterly MPC (%)'
				            'Annual MPC (%)'
				            'Beta (Annualized)'
				};
			new_entries = {	round(values.mpcs(5).avg_quarterly * 100, 1)
		                    round(values.mpcs(5).avg_annual * 100, 1)
		                    round(values.beta_annualized, 3) 
				};
			out = tables.TableGen.append_to_table(out,...
				new_entries, new_labels);
		end
	end
end