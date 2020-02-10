clear

addpath('/home/brian/Documents/GitHub/Discrete_HA/')
addpath('/home/brian/Documents/GitHub/Discrete_HA/code')

load('/home/brian/Documents/GitHub/Discrete_HA/output/variables1.mat')
load('/home/brian/Documents/GitHub/Discrete_HA/input/agrid.mat')
load('/home/brian/Documents/GitHub/Discrete_HA/input/yPdist.mat')

params = Sparams;
mpc_plotter = statistics.MPCPlotter(params, agrid, yPdist, results);


%% Plot MPCs
yP_indices = [3, 8];
mpc_plotter.plot_mpcs(yP_indices)

%% Plot wealth histogram
mpc_plotter.plot_asset_dist(200, 10)