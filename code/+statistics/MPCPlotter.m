classdef MPCPlotter < handle

	properties
		p;
		grids;
		dims;
		mpcs;
		agrid;
		pmf_a;

		lims = [0, 10, -inf, inf];
		fig;
		ax;
		tile;

		mpcs_plotted = false;

		fontsize = 12;
		show_grid = 'off';
	end

	methods
		function obj = MPCPlotter(params, agrid, results)
			obj.p = params;
			obj.dims = [numel(agrid), obj.p.nyP, obj.p.nyF, obj.p.nb];

			obj.mpcs = reshape(...
				results.direct.mpcs(5).mpcs_1_t{1}, obj.dims);

			obj.agrid = agrid;
			obj.pmf_a = results.direct.agrid_dist;
		end

		function [ax_main, ax_window] = create_mpcs_plot(...
			obj, yP_indices, zoomed_window, shock_size)

			obj.fig = figure();
			ax_main = axes(obj.fig);
			legend('hide');
			xlabel("Wealth (ratio to mean annual income)")
			ylabel(sprintf("MPC out of %g", shock_size))
			ax_main = obj.apply_formatting(ax_main);

			if zoomed_window
				ax_window = add_window(obj.fig, ax_main);
				ax_window = obj.apply_formatting(ax_window);
			end

			for ii = 1:numel(yP_indices)
			    iyP = yP_indices(ii);
			    
			    obj.plot_mpc_function(ax_main, iyP);
			    if zoomed_window
			        obj.plot_mpc_function(ax_window, iyP);
			    end
			end

			yP_labels = {'Low yP', 'High yP'};
			legend(ax_main, yP_labels);
			legend(ax_main, 'show');
			box(ax_main.Legend, 'off');
		end

		function plot_mpc_function(obj, parent_obj, iyP, iyF, ib)
			if nargin == 3
				iyF = 1;
				ib = 1;
			elseif nargin == 4
				ib = 1;
			end

			selected_mpcs = obj.mpcs(:,iyP,iyF,ib);

			hold(parent_obj, 'on');
			plot(parent_obj, obj.agrid, selected_mpcs);
			hold(parent_obj, 'off');
		end

		function ax = add_median_wealth(obj, ax, median_wealth)
			hold(ax, 'on')
			wline = [0, 1];
			wvals = [median_wealth, median_wealth];
			plot(ax, wvals, wline, '--b');
		    ax.Legend.String{end} = 'Median wealth';
		    hold(ax, 'off')
		end

		function ax = apply_formatting(obj, ax)
			grid(ax, obj.show_grid);
			set(ax, 'FontSize', obj.fontsize);
		end
	end
end

function ax = add_window(parent_obj, main_ax)
	main_ax_corner = main_ax.Position(1:2);
	main_ax_width = main_ax.Position(3);
	main_ax_height = main_ax.Position(4);

	window_corner(1) = main_ax_corner(1) + 0.25 * main_ax_width;
	window_corner(2) = main_ax_corner(2) + 0.4 * main_ax_height;
	window_width = 0.5 * main_ax_width;
	window_height = 0.3 * main_ax_height;

	window_position = [window_corner, window_width, window_height];
	ax = axes('Position', window_position);
	box(ax, 'on');
end
