classdef TableGenerator
	properties (SetAccess = private)
		output_table;
	end

	properties
		decomp_baseline;
		decomp_incrisk_alt;
	end

	methods
		function output_table = create(obj, params, results, freq)
			output_table = table();
			this_freq = find([params.freq]==freq);
			if isempty(this_freq)
			    return;
			end

			shock_labels = params(1).shocks_labels;

			mpcs_present = false;
			mpcs_news_present = false;
			mpcs_loan_loss_present = false;
			decomp_norisk_present = false;
			decomp_RA_present = false;
			for ip = this_freq
				if params(ip).MPCs
					mpcs_present = true;
				end

				if params(ip).MPCs_news
					mpcs_news_present = true;
				end

				if params(ip).MPCs_loan_and_loss
					mpcs_loan_loss_present = true;
				end

				if results(ip).decomp_norisk(1).completed
					decomp_norisk_present = true;
				end

				if results(ip).decomp_RA.completed
					decomp_RA_present = true;
				end
			end

			for ip = this_freq
				p = params(ip);
				result_structure = results(ip);

				new_column = intro_stats_table(result_structure, p);

				temp = wealth_constraints_table(result_structure);
				new_column = [new_column; temp];

				temp = wealth_distribution_table(result_structure);
				new_column = [new_column; temp];

				ishock = 1;
				while mpcs_present && (ishock < 7)
					shock_size = shock_labels{ishock};
					temp = mpcs_table(result_structure, p, ishock, shock_size);
					new_column = [new_column; temp];
					ishock = ishock + 1;
				end

				if 	~isempty(decomp_norisk_present)
					shock_size = shock_labels{5};

					decomp_structure = results(ip).decomp_norisk;
					temp = decomp_table(decomp_structure, shock_size);
					new_column = [new_column; temp];
				end

				if ~isempty(decomp_RA_present)
					shock_size = shock_labels{5};

					decomp_structure = results(ip).decomp_RA;
					temp = repagent_decomp_table(decomp_structure, shock_size);
					new_column = [new_column; temp];
				end

				if ~isempty(obj.decomp_baseline)
					shock_size = shock_labels{5};

					decomp_structure = obj.decomp_baseline(ip);
					temp = baseline_decomp_table(decomp_structure, shock_size);
					new_column = [new_column; temp];
				end

				if ~isempty(obj.decomp_incrisk_alt)
					shock_size = shock_labels{5};

					decomp_structure = obj.decomp_incrisk_alt(ip);
					temp = alt_decomp_table(decomp_structure, shock_size);
					new_column = [new_column; temp];
				end

				if mpcs_news_present
					temp = mpcs_news_table(result_structure, p, 'MEAN',...
						shock_labels);
					new_column = [new_column; temp];

					temp = mpcs_news_table(result_structure, p, 'MEDIAN',...
						shock_labels);
					new_column = [new_column; temp];
				end

				if mpcs_loan_loss_present
					temp = specialty_mpc_tables(result_structure, p,...
						shock_labels);
					new_column = [new_column; temp];
				end

				column_label = sprintf('Specification%d', p.index);
				new_column.Properties.VariableNames = {column_label};
				output_table = [output_table, new_column];
			end
		end
	end
end

function out = intro_stats_table(values, p)
	out = table({p.name},...
		'VariableNames', {'results'},...
		'RowNames', {'Model'});

	new_labels = {	'Beta (Annualized)'
		            'Mean gross annual income'
		            'Stdev log annual gross income'
		            'Stdev log annual net income'
		};
	new_entries = {	values.direct.beta_annualized
                    values.direct.mean_grossy_A
                    values.direct.stdev_loggrossy_A
                    values.direct.stdev_lognety_A    
		};
	out = append_to_table(out, new_entries, new_labels);
end

function out = wealth_constraints_table(values)
	out = new_table_with_header('WEALTH CONSTRAINTS');

	% Mean assets and saving
	new_labels = {	'Mean assets'
					'Median assets'
					'Mean saving'
		};
	new_entries = {	values.direct.mean_a
					values.direct.median_a
					values.direct.mean_s
		};
	out = append_to_table(out, new_entries, new_labels);

	% Fraction with saving = 0
	new_labels = {'Fraction with s == 0'};
	new_entries = values.direct.s0;
	out = append_to_table(out, new_entries, new_labels);

	% Fraction with assets or cash<= some value
	new_labels = {	'Fraction with a == 0'
		            'Fraction with a <= 0.5% mean ann gross lab inc'
		            'Fraction with a <= 1% mean ann gross lab inc'
		            'Fraction with a <= 2% mean ann gross lab inc'
		            'Fraction with a <= 5% mean ann gross lab inc'
		            'Fraction with a <= 10% mean ann gross lab inc'
		            'Fraction with a <= 15% mean ann gross lab inc'
		            'Fraction with a <= $1000'
		            'Fraction with a <= $5000'
		            'Fraction with a <= $10000'
		            'Fraction with a <= 1/6 own quarterly income'
		            'Fraction with a <= 1/12 own quarterly income'
		            'Fraction with x <= 1/6 own quarterly income'
		            'Fraction with x <= 1/12 own quarterly income'
		};
	temp = num2cell(values.direct.constrained(:));
	new_entries = {	values.direct.wealth_lt_1000
	                values.direct.wealth_lt_5000
	                values.direct.wealth_lt_10000
	                values.direct.a_sixth_sim
	                values.direct.a_twelfth_sim
	                values.direct.x_sixth_sim
                	values.direct.x_twelfth_sim
		};
	new_entries = [temp; new_entries];
	out = append_to_table(out, new_entries, new_labels);
end

function out = wealth_distribution_table(values)
	out = new_table_with_header('WEALTH DISTRIBUTION');

	% Percentiles
	new_labels = {	'Wealth, 10th percentile'
		            'Wealth, 25th percentile'
		            'Wealth, 50th percentile'
		            'Wealth, 75th percentile'
		            'Wealth, 90th percentile'
		            'Wealth, 95th percentile'
		            'Wealth, 99th percentile'
		            'Wealth, 99.9th percentile'
		};
	new_entries = num2cell(values.direct.wpercentiles(:));
	out = append_to_table(out, new_entries, new_labels);

	% Other stats
	new_labels = {	'Wealth, top 10% share'
					'Wealth, top 1% share'
					'Gini coefficient'
                    'Rank-rank correlation, assets and pref het'
		};
	new_entries = {	values.direct.top10share
					values.direct.top1share
					values.direct.wealthgini
                    values.direct.assets_z_rank_corr
		};
	out = append_to_table(out, new_entries, new_labels);
end

function out = mpcs_table(values, p, ishock, shock_size)
	header = sprintf(...
		'MPCs OUT OF FRACTION MEAN ANN INC (SHOCK %s)', shock_size);
	out = new_table_with_header(header);

	new_labels = {	sprintf('PERIOD 1 Median(MPC), shock = %s', shock_size)
					sprintf('PERIOD 1 E[MPC], shock = %s', shock_size)
					sprintf('PERIOD 2 E[MPC], shock = %s', shock_size)
					sprintf('PERIOD 3 E[MPC], shock = %s', shock_size)
					sprintf('PERIOD 4 E[MPC], shock = %s', shock_size)
					sprintf('FOUR PERIOD E[MPC], shock = %s', shock_size)
		};

	if p.MPCs
		new_entries = {	values.direct.mpcs(ishock).median(1,1);
						values.direct.mpcs(ishock).avg_s_t(1,1)
						values.direct.mpcs(ishock).avg_s_t(1,2)
						values.direct.mpcs(ishock).avg_s_t(1,3)
						values.direct.mpcs(ishock).avg_s_t(1,4)
						values.direct.mpcs(ishock).avg_1_1to4
			};
	else
		new_entries = num2cell(NaN(6,1));
	end
	out = append_to_table(out, new_entries, new_labels);
end

function out = mpcs_news_table(values, p, statistic, shock_labels)
	if strcmp(statistic, 'MEAN')
		stat_label = 'E[MPC]';
	elseif strcmp(statistic, 'MEDIAN')
		stat_label = 'Median(MPC)';
	else
		error('invalid statistic');
	end

	header_name = strcat(statistic, ' MPCs OUT OF NEWS');
	out = new_table_with_header(header_name);

	new_labels = {};
	new_entries = {};
	ilabel = 1;
	for period = [2, 5]
		for ishock = 1:6
			shock_size = shock_labels{ishock};
			shock_label = sprintf(', shock = %s in %d period(s)', shock_size, period-1);

			new_labels{ilabel} = convertStringsToChars(...
				strcat("ONE PERIOD ", stat_label, shock_label));

			if ~p.MPCs_news
				new_entries{ilabel} = NaN;
			elseif strcmp(statistic, 'MEAN')
				new_entries{ilabel} = values.direct.mpcs(ishock).avg_s_t(period,1);
			elseif strcmp(statistic, 'MEDIAN')
				new_entries{ilabel} = values.direct.mpcs(ishock).median(period,1);
			end
			ilabel = ilabel + 1;
		end
	end

	new_labels = reshape(new_labels, [], 1);
	new_entries = reshape(new_entries, [], 1);

	out = append_to_table(out, new_entries, new_labels);
end

function out = specialty_mpc_tables(values, p, shock_labels)
	% Loss in two years
	loss_size = shock_labels{1};
	header_name = sprintf('MPCs OUT OF LOSS IN 8 PERIODS (SHOCK %s)', loss_size);
	loss_table = new_table_with_header(header_name);

	new_entries = {	values.direct.loss_in_2_years.avg
                    values.direct.loss_in_2_years.mpc_condl
                    values.direct.loss_in_2_years.mpc_neg
                    values.direct.loss_in_2_years.mpc0
                    values.direct.loss_in_2_years.mpc_pos
                    values.direct.loss_in_2_years.median
        };
    shock_label = sprintf(', shock = %s in 8 period(s)', loss_size);
    new_labels = {	strcat('ONE PERIOD E[MPC]', shock_label)
    				strcat('ONE PERIOD E[MPC|MPC>0]', shock_label)
    				strcat('ONE PERIOD P(MPC<0)', shock_label)
    				strcat('ONE PERIOD P(MPC=0)', shock_label)
    				strcat('ONE PERIOD P(MPC>0)', shock_label)
    				strcat('ONE PERIOD Median(MPC)', shock_label)
    	};
    loss_table = append_to_table(loss_table, new_entries, new_labels);

	% One-year loan
	loan_size = shock_labels{6};
	header_name = sprintf('MPCs OUT OF 4-PERIOD LOAN OF %s', loan_size);
	loan_table = new_table_with_header(header_name);

	new_entries = {	values.direct.loan.avg
                    values.direct.loan.mpc_condl
                    values.direct.loan.mpc_neg
                    values.direct.loan.mpc0
                    values.direct.loan.mpc_pos
                    values.direct.loan.median
        };
    shock_label = sprintf(', loan of %s', loan_size);
    new_labels = {	strcat('ONE PERIOD E[MPC]', shock_label)
    				strcat('ONE PERIOD E[MPC|MPC>0]', shock_label)
    				strcat('ONE PERIOD P(MPC<0)', shock_label)
    				strcat('ONE PERIOD P(MPC=0)', shock_label)
    				strcat('ONE PERIOD P(MPC>0)', shock_label)
    				strcat('ONE PERIOD Median(MPC)', shock_label)
    	};
    loan_table = append_to_table(loan_table, new_entries, new_labels);

    out = [loss_table; loan_table];
end

function out = decomp_table(decomp, shock_size)
	header_name = sprintf(...
		'DECOMP OF ONE PERIOD E[MPC], SHOCK OF %s', shock_size);
	out = new_table_with_header(header_name);

	new_labels = {	'Decomp of E[MPC] around a=0, RA MPC'
		            'Decomp of E[MPC] around a=0, HtM Effect'
		            'Decomp of E[MPC] around a=0, Non-HtM, constraint'
		            'Decomp of E[MPC] around a=0, Non-HtM, inc risk'
		            'Decomp of E[MPC] around a=0.01, RA MPC'
		            'Decomp of E[MPC] around a=0.01, HtM Effect'
		            'Decomp of E[MPC] around a=0.01, Non-HtM, constraint'
		            'Decomp of E[MPC] around a=0.01, Non-HtM, inc risk'
		            'Decomp of E[MPC] around a=0.05, RA MPC'
		            'Decomp of E[MPC] around a=0.05, HtM Effect'
		            'Decomp of E[MPC] around a=0.05, Non-HtM, constraint'
		            'Decomp of E[MPC] around a=0.05, Non-HtM, inc risk'
		};

	temp = [	[decomp.term1]
				[decomp.term2]
				[decomp.term3]
				[decomp.term4]
		];
	new_entries = num2cell(temp(:));

	out = append_to_table(out, new_entries, new_labels);
end

function out = repagent_decomp_table(decomp, shock_size)
	header_name = sprintf(...
		'DECOMP OF ONE PERIOD E[MPC] - MPC_RA, SHOCK OF %s', shock_size);
	out = new_table_with_header(header_name);

	new_labels = {	'E[MPC] - MPC_RA'
		            'Decomp of E[MPC] - MPC_RA, effect of MPC fcn'
		            'Decomp of E[MPC] - MPC_RA, effect of distr'
		            'Decomp of E[MPC] - MPC_RA, interaction'
		};
	new_entries = {	decomp.Em1_less_mRA
                    decomp.term1
                    decomp.term2
                    decomp.term3
		};

	out = append_to_table(out, new_entries, new_labels);
end

function out = baseline_decomp_table(decomp, shock_size)
	header_name = sprintf(...
		'DECOMP OF EM1-EM0, SHOCK OF %s', shock_size);
	out = new_table_with_header(header_name);

	new_labels = {	'Em1 - Em0'
		            'Decomp of Em1-Em0, effect of MPC fcn'
		            'Decomp of Em1-Em0, effect of distr'
		            'Decomp of Em1-Em0, interaction'
		            'Effect of MPC fcn, level'
		            'Effect of MPC fcn, shape'
		            'Effect of distr around 0, HtM households'
		            'Effect of distr around 0, non-HtM households'
		            'Effect of distr around 0.01, HtM households'
		            'Effect of distr around 0.01, non-HtM households'
		            'Effect of distr around 0.05, HtM households'
		            'Effect of distr around 0.05, non-HtM households'
		};

	new_entries = {	decomp.Em1_less_Em0
                    decomp.term1                       
                    decomp.term2
                    decomp.term3
                    decomp.term1a
                    decomp.term1b
                    decomp.term2a(1)   
                    decomp.term2b(1)
                    decomp.term2a(2)   
                    decomp.term2b(2)
                    decomp.term2a(3)   
                    decomp.term2b(3)
        };
	out = append_to_table(out, new_entries, new_labels);
end

function out = alt_decomp_table(decomp, shock_size)
	header_name = sprintf(...
		'ALT DECOMP OF ONE PERIOD E[MPC], SHOCK OF %s', shock_size);
	out = new_table_with_header(header_name);

	new_labels = {	'Alt decomp of E[MPC] around a=0, RA MPC'
		            'Alt decomp of E[MPC] around a=0, HtM Effect'
		            'Alt decomp of E[MPC] around a=0, Non-HtM, constraint'
		            'Alt decomp of E[MPC] around a=0, Non-HtM, inc risk'
		            'Alt decomp of E[MPC] around a=0, Interaction'
		            'Alt decomp of E[MPC] around a=0.01, RA MPC'
		            'Alt decomp of E[MPC] around a=0.01, HtM Effect'
		            'Alt decomp of E[MPC] around a=0.01, Non-HtM, constraint'
		            'Alt decomp of E[MPC] around a=0.01, Non-HtM, inc risk'
		            'Alt decomp of E[MPC] around a=0.01, Interaction'
		            'Alt decomp of E[MPC] around a=0.05, RA MPC'
		            'Alt decomp of E[MPC] around a=0.05, HtM Effect'
		            'Alt decomp of E[MPC] around a=0.05, Non-HtM, constraint'
		            'Alt decomp of E[MPC] around a=0.05, Non-HtM, inc risk'
		            'Alt decomp of E[MPC] around a=0.05, Interaction'
		};

	temp = [	[decomp.term1]
				[decomp.term2]
				[decomp.term3]
				[decomp.term4]
				[decomp.term5]
		];
	new_entries = num2cell(temp(:));

	out = append_to_table(out, new_entries, new_labels);
end

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