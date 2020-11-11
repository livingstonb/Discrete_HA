classdef OtherPanels
	methods (Static)
		function out = intro_panel(values, p, group)
			if isempty(p.label)
				param_label = p.name;
			else
				param_label = p.label;
			end

			if nargin < 3
				group = '';
			end

			out = table({param_label},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			if strcmp(group, 'Q2')
				new_labels = {	'Spacing'
								'Switching probability'
					};
				new_entries = {	p.betawidth
								p.prob_zswitch
					};
				out = tables.TableGen.append_to_table(out,...
					new_entries, new_labels);
			elseif ismember(group, {'Q3', 'Q4'})
				new_labels = {	'Risk aversion'
								'IES'
					};

				if isempty(p.other)
					new_entries = {	p.risk_aver
									1 / p.invies
						};
					new_entries = aux.cellround(new_entries, 3);
				else
					new_entries = { p.other{1}
									p.other{2}
						};
				end

				out = tables.TableGen.append_to_table(out,...
					new_entries, new_labels);
			elseif strcmp(group, 'Q5')
				new_labels = {'Temptation Parameter'};

				if isempty(p.other)
					new_entries = {p.temptation};
				else
					new_entries = { p.other{1} };
				end
				out = tables.TableGen.append_to_table(out,...
					new_entries, new_labels);
			elseif strcmp(group, 'Q6')
				new_labels = {'Annualized returns (%)'};

				if isempty(p.other)
					new_entries = { (p.R^(p.freq) - 1) * 100 };
				else
					new_entries = { p.other{1} };
				end
				out = tables.TableGen.append_to_table(out,...
					new_entries, new_labels);
			elseif strcmp(group, 'Q7')
				new_labels = {'Description'};
				new_entries = { p.other{1} };
				out = tables.TableGen.append_to_table(out,...
					new_entries, new_labels);
			end

			new_labels = {	'Quarterly MPC (%)'
				            'Annual MPC (%)'
				            'Beta (Annualized)'
				};
			new_entries = {	round(values.mpcs(5).quarterly.value * 100, 1)
		                    round(values.mpcs(5).annual.value * 100, 1)
		                    round(values.beta_annualized, 3) 
				};
			out = tables.TableGen.append_to_table(out,...
				new_entries, new_labels);
		end
	end
end