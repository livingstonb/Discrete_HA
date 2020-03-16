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

				% shock = stats_ip.mpcs(5).shock.value;
				% obj.mpc_comparison(stats_ip, shock);
				% obj.mpc_comparison_pct(stats_ip, shock);

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

			new_entries = {
				stats.mpcs(5).quarterly
				stats.mpcs(5).annual
				stats.beta_A
				stats.mean_gross_y_annual
				stats.std_log_gross_y_annual
				stats.std_log_net_y_annual
			};

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
				stats.a_lt_ysixth
				stats.a_lt_ytwelfth
				stats.w_top10share
				stats.w_top1share
				stats.wgini
			};

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
			obj.update_current_column(out, new_entries);
		end

		function decomp_norisk_table(obj, p, stats)
			panel_name = 'Decomps of E[MPC] wrt RA and no inc risk, $500 shock';
			out = obj.new_table_with_header(panel_name);

			tmp = stats.mpcs(5).quarterly;
			tmp.label = 'Quarterly MPC (%)'
			new_entries = {tmp, stats.decomp_norisk.term1_pct};
			obj.update_current_column(out, new_entries);

			for ithresh = 1:numel(p.abars)
				threshold = p.abars(ithresh);
				panel_name = sprintf('For HtM threshold #%d', ithresh);
				out = obj.new_table_with_header(panel_name);

				new_entries = {
					stats.decomp_norisk.term2(ithresh);
					stats.decomp_norisk.term3(ithresh);
					stats.decomp_norisk.term4(ithresh);
				};

				obj.update_current_column(out, new_entries);
			end
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
			}

			obj.update_current_column(out, new_entries);
		end
	end
end