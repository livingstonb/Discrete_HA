clear
close all

addpath('/home/brian/Documents/GitHub/Discrete_HA/')
addpath('/home/brian/Documents/GitHub/Discrete_HA/code')

load('/home/brian/Documents/GitHub/Discrete_HA/output/variables1.mat')
load('/home/brian/Documents/GitHub/Discrete_HA/input/agrid.mat')
load('/home/brian/Documents/GitHub/Discrete_HA/input/yPdist.mat')

params = Sparams;

fontsize = 16;

%% MPCs Function
mpc_plotter = statistics.MPCPlotter(params, agrid, yPdist, results);
mpc_plotter.fontsize = fontsize;
mpc_plotter.show_grid = 'on';

yP_indices = [3, 8];
zoomed_window = true;
shock_size = 0.01;
[ax_main, ax_window] = mpc_plotter.create_mpcs_plot(...
			yP_indices, zoomed_window, shock_size);
ax_main.XLim = [0, 10];

window_max_x = 0.3;
ax_window.YLim = ax_main.YLim;
ax_window.XLim = [0, window_max_x];
xticks(ax_window, [0:0.1:window_max_x])
yticks(ax_window, [0:0.1:0.3])
set(ax_window, 'FontSize', fontsize-2)
        
%% Plot wealth histogram
mpc_plotter = statistics.MPCPlotter(params, agrid, yPdist, results);
[ax, wealth_hist] = mpc_plotter.create_wealth_histogram(200, 10);

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