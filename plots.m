clear
close all

addpath('/home/brian/Documents/GitHub/Discrete_HA/')
addpath('/home/brian/Documents/GitHub/Discrete_HA/code')

outdir = '/home/brian/Documents/GitHub/Discrete_HA/output';

load('/home/brian/Documents/GitHub/Discrete_HA/output/variables.mat')
load('/home/brian/Documents/GitHub/Discrete_HA/input/agrid_nolin.mat')
load('/home/brian/Documents/GitHub/Discrete_HA/input/yPdist.mat')

params = Sparams;

fontsize = 12;

%% Policy Function

%% MPCs Function
% mpc_plotter = statistics.MPCPlotter(params, agrid, results);
% mpc_plotter.fontsize = fontsize;
% mpc_plotter.show_grid = 'on';
% 
% yP_indices = [3, 8];
% zoomed_window = true;
% shock_size = 0.01;
% [ax_main, ax_window] = mpc_plotter.create_mpcs_plot(...
% 			yP_indices, zoomed_window, shock_size);
% ylim_main = ax_main.YLim;
%         
% imedian = find(params.percentiles == 50);
% median_wealth = results.direct.wpercentiles(imedian);
% ax_main = mpc_plotter.add_median_wealth(ax_main, median_wealth);
% 
% ax_main.XLim = [0, 5];
% ax_main.YLim = ylim_main;
% 
% window_max_x = 0.3;
% ax_window.YLim = ax_main.YLim;
% ax_window.XLim = [0, window_max_x];
% xticks(ax_window, [0:0.1:window_max_x])
% yticks(ax_window, [0:0.1:0.3])
% set(ax_window, 'FontSize', fontsize-2)
% ax_window.YTick = ax_main.YTick(1:2:end);
% 
% figpath = fullfile(outdir, 'mpc_function.jpg');
% saveas(gcf, figpath)

%% Wealth pmf
plot(agrid, results.direct.agrid_dist)
xlim([0 0.2])
        
%% Wealth histogram
nbins = 100;
amax = {0.25};
amax_visible = 0.2;

iyP = 1:11;
nyP = numel(iyP);

pmf_a = results.direct.adist(:,iyP,:,:);
pmf_a = pmf_a(:) / sum(pmf_a(:));
pmf_a = reshape(pmf_a, [], nyP);
pmf_a = sum(pmf_a, 2);

            
wealth_plotter = statistics.WealthPlotter(params, agrid, pmf_a);
[ax, wealth_hist] = wealth_plotter.create_histogram(nbins, amax{:});
title("Wealth condl on low yP, truncated above")
ax.XLim = [0, amax_visible];
ax.YLim = [0, max(wealth_hist.Values(1:end-1))];

figpath = fullfile(outdir, 'wealth_condl_on_low_yP.jpg');
saveas(gcf, figpath)


% %% MPCs Function
% mpc_plotter = statistics.MPCPlotter(params, agrid, yPdist, results);
% mpc_plotter.plot_zoomed_window = true;
% mpc_plotter.font_size = 16;
% 
% fig = figure();
% 
% ax_main = statistics.MPCPlotter.gen_formatted_axis(fig, 0.01);
% ax_main.FontSize = 16;
% 
% if plot_zoomed_window
%     ax_window = statistics.MPCPlotter.add_window(fig, ax_main);
%     ax_window.FontSize = 18;
% end
% 
% yP_indices = [3, 8];
% yP_labels = {'Low', 'High'};
% 
% for ii = 1:numel(yP_indices)
%     iyP = yP_indices(ii);
%     label = yP_labels{ii};
%     
%     mpc_plotter.plot_mpc_function(ax_main, iyP);
%     mpc_plotter.set_new_legend_entry(ax_main, label);
%     if plot_zoomed_window
%         mpc_plotter.plot_mpc_function(ax_window, iyP);
%     end
% end
% 
% statistics.MPCPlotter.format_picture_in_picture(ax_main, ax_window)

% %% Plot MPCs
% yP_indices = [3, 8];
% mpc_plotter.plot_mpcs(yP_indices)
% 
% %% Plot wealth histogram
% mpc_plotter.plot_asset_dist(200, 10)