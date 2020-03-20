classdef StatsTable < tables.BaseTable
	methods
		function output_table = create(obj, params, stats)
			obj.output = table();
			for ip = 1:obj.n_cols
				p_ip = params(ip);
				stats_ip = stats{ip};

				obj.current_column = table();
				obj.intro_stats_table(p_ip, stats_ip);
				obj.wealth_stats_table(stats_ip);
				obj.mpc_size_table(stats_ip);
				obj.mpc_sign_table(stats_ip);

				obj.mpc_comparison(ip);
				obj.mpc_comparison_pct(ip);

				obj.percentiles_tables(stats_ip);

    %             obj.other_params_table(stats_ip);
    			obj.decomp_ra(stats_ip);
    			obj.decomp_norisk_table(p_ip, stats_ip);
                obj.other_stats_table(stats_ip);

                obj.add_column(ip)
			end

			output_table = obj.output;
		end

		function intro_stats_table(obj, p, stats)
			out = table({p.label},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			descr = struct();
			if isempty(p.other)
				descr.value = '';
			else
				descr.value = p.other{1};
			end
			descr.label = 'Description';

			new_entries1 = {
				stats.params.freq
				descr
			};

			new_entries2 = {
				stats.mpcs(5).quarterly
				stats.mpcs(5).annual
			};

			new_entries3 = {
				stats.beta_A
				stats.mean_gross_y_annual
				stats.std_log_gross_y_annual
				stats.std_log_net_y_annual
			};


			new_entries2 = obj.sround_mult(new_entries2, 1);
			new_entries3 = obj.sround_mult(new_entries3, 3);
			new_entries = [new_entries1; new_entries2; new_entries3];

			obj.update_current_column(out, new_entries);
		end

		function wealth_stats_table(obj, stats)
			panel_name = 'Wealth Statistics';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.mean_a
				stats.median_a
				stats.sav0
				stats.constrained_pct{1}
				stats.constrained_pct{2}
				stats.constrained_pct{3}
				stats.constrained_pct{4}
				stats.constrained_pct{5}
				stats.constrained_pct{6}
				stats.constrained_pct{7}
				stats.a_lt_ysixth
				stats.a_lt_ytwelfth
				stats.w_top10share
				stats.w_top1share
				stats.wgini
			};

			new_entries = obj.sround_mult(new_entries, 3);

			obj.update_current_column(out, new_entries);
		end

		function mpc_size_table(obj, stats)
			panel_name = 'MPC size effects';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(4).quarterly
				stats.mpcs(6).quarterly
				stats.mpcs(4).annual
				stats.mpcs(6).annual
			};

			new_entries = obj.sround_mult(new_entries, 1);

			obj.update_current_column(out, new_entries);
		end

		function mpc_sign_table(obj, stats)
			panel_name = 'MPC sign effects';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(1).quarterly
				stats.mpcs(2).quarterly
				stats.mpcs(3).quarterly
				stats.mpcs(1).annual
				stats.mpcs(2).annual
				stats.mpcs(3).annual
			};

			new_entries = obj.sround_mult(new_entries, 1);

			obj.update_current_column(out, new_entries);
		end

		function percentiles_tables(obj, stats)
			panel_name = 'Wealth Distribution';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.wpercentiles{1}
				stats.wpercentiles{2}
				stats.wpercentiles{3}
				stats.wpercentiles{5}
				stats.wpercentiles{7}
				stats.wpercentiles{8}
			};

			new_entries = obj.sround_mult(new_entries, 3);

			obj.update_current_column(out, new_entries);
		end

		function decomp_ra(obj, stats)
			panel_name = 'Decomposition of E[MPC] - MPC_RA';
			out = obj.new_table_with_header(panel_name);

			new_entries = {	stats.decomp_RA.Em1_less_mRA
							stats.decomp_RA.term1
		                    stats.decomp_RA.term2
		                    stats.decomp_RA.term3
				};

			new_entries = obj.sround_mult(new_entries, 3);

			obj.update_current_column(out, new_entries);
		end

		function decomp_norisk_table(obj, p, stats)
			panel_name = 'Decomps of E[MPC] wrt RA and no inc risk, $500 shock';
			out = obj.new_table_with_header(panel_name);

			if p.freq == 1
				tmp = stats.mpcs(5).annual;
			else
				tmp = stats.mpcs(5).quarterly;
			end
			tmp.label = 'MPC (%)';
			new_entries = {
				tmp
				stats.decomp_norisk.term1_pct
				};
			new_entries = obj.sround_mult(new_entries, 3);
			obj.update_current_column(out, new_entries);

			for ithresh = 1:numel(p.abars)
				threshold = p.abars(ithresh);
				panel_name = sprintf('For HtM threshold #%d', ithresh);
				out = obj.new_table_with_header(panel_name);

				new_entries = {
					stats.decomp_norisk.term2(ithresh)
					stats.decomp_norisk.term3(ithresh)
					stats.decomp_norisk.term4(ithresh)
				};

				new_entries = obj.sround_mult(new_entries, 3);

				obj.update_current_column(out, new_entries);
			end
		end

		function mpc_comparison(obj, ip)
			if isempty(obj.decomp_baseline)
				return
			end

			decomp_ip = obj.decomp_baseline(ip);
			panel_name = decomp_ip.description.value;
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				decomp_ip.Em1_less_Em0
				decomp_ip.term1
				decomp_ip.term1a
				decomp_ip.term1b
				decomp_ip.term2
				decomp_ip.term2a(3)
				decomp_ip.term2b(3)
				decomp_ip.term3
			};

			new_entries = obj.sround_mult(new_entries, 3);

			obj.update_current_column(out, new_entries);
		end

		function mpc_comparison_pct(obj, ip)
			if isempty(obj.decomp_baseline)
				return
			end

			decomp_ip = obj.decomp_baseline(ip);
			panel_name = 'Decomposition in terms of % of E[MPC] - E[MPC_b]';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				decomp_ip.term1_pct
				decomp_ip.term1a_pct
				decomp_ip.term1b_pct
				decomp_ip.term2_pct
				decomp_ip.term2a_pct(3)
				decomp_ip.term2b_pct(3)
				decomp_ip.term3_pct
			};

			new_entries = obj.sround_mult(new_entries, 1);

			obj.update_current_column(out, new_entries);
		end

		function other_stats_table(obj, stats)
			panel_name = 'Other Statistics';
			out = obj.new_table_with_header(panel_name);

			new_entries = {
				stats.constrained_dollars{1}
				stats.constrained_dollars{2}
				stats.constrained_dollars{3}
				stats.constrained_dollars{4}
				stats.constrained_dollars{5}
			};

			new_entries = obj.sround_mult(new_entries, 3);

			obj.update_current_column(out, new_entries);
		end
	end
end