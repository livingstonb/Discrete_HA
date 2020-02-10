classdef MPCPlotter < handle

	properties
		p;
		grids;
		dims;
		yPdist;
		mpcs;
		agrid;
		pmf_a;

		lims = [0, 10, -inf, inf];
		fig;
		ax;
		tile;

		mpcs_plotted = false;
	end

	methods
		function obj = MPCPlotter(params, agrid, yPdist, results)
			obj.p = params;
			obj.dims = [obj.p.nx_DST, obj.p.nyP, obj.p.nyF, obj.p.nb];
			obj.yPdist = yPdist;

			obj.mpcs = reshape(...
				results.direct.mpcs(5).mpcs_1_t{1}, obj.dims);

			obj.agrid = agrid;
			obj.pmf_a = results.direct.agrid_dist;

			obj.fig = figure();
			
			obj.tile = tiledlayout(2, 1);
			% obj.ax = axes('Parent', obj.fig);
			% hold(obj.ax, 'on');
		end

		function plot_mpcs(obj, yP_indices, iyF, ib)
			if nargin == 2
				iyF = 1;
				ib = 1;
			elseif nargin == 2
				ib = 1;
			end

			ax = nexttile(obj.tile);
			hold(ax, 'on');
			for iyP = yP_indices
				selected_mpcs = obj.mpcs(:,iyP,iyF,ib);
				plot(ax, obj.agrid, selected_mpcs);
				% hold on
			end
			hold(ax, 'off');

			legend("Low persistent income", "High persistent income")
			axis(obj.lims);
			xlabel("Wealth (ratio to mean annual income)")
			ylabel("One-period MPC out of 0.01")

			obj.mpcs_plotted = true;
		end

		function plot_asset_dist(obj, nbins, amax)
			% if obj.mpcs_plotted
			% 	obj.ax.OuterPosition = [0 0.5 1 1];
			% 	yyaxis(obj.ax, 'right');
			% 	obj.ax.OuterPosition = [0 0 1 1];
			% end
			ax = nexttile(obj.tile);
			[edges, counts] = smoothed_histogram(obj.agrid, obj.pmf_a, nbins, amax);
			asset_hist = histogram('BinEdges', edges, 'BinCounts', counts);

			bar_color = [0,0,0] + 0.5;
			asset_hist.FaceColor = bar_color;
			% asset_hist.FaceAlpha = 0.2;
			asset_hist.EdgeColor = bar_color;
			xlabel("Wealth (ratio to mean annual income)")
			ylabel("Probability density")

			% bfig = bar(ax, bins, vals, 'hist');

			% bar_color = [0,0,0] + 0.5;
			% bfig.FaceColor = bar_color;
			% bfig.FaceAlpha = 0.2;
			% bfig.EdgeColor = bar_color;

			% yyaxis(obj.ax, 'left');
		end
	end
end

function [bins, vals] = smoothed_histogram(agrid, pmf, nbins, amax)
	% dims = size(pmf);
	% dims_for_rep = dims;
	% dims_for_rep(1) = 1;

	% agrid_rep

	% n_asset = numel(agrid);
	% pmf_reshaped = reshape(pmf, n_asset, []);
	% pmf_marginal = sum(pmf, 2);
	% tmp_sorted = sortrows([agrid(:), pmf(:)]);
	% tmp_sorted = sortrows([agrid, pmf_marginal]);

	% [a_vals, inds_unique] = unique(tmp_sorted(:,1), 'last');
	% a_cdf = cumsum(tmp_sorted(:,2));
	% a_cdf = a_cdf(inds_unique);

	% cdf_interp = griddedInterpolant(a_vals, a_cdf, 'linear');

	% a_cdf = cumsum(pmf_marginal);
	a_cdf = cumsum(pmf);
	cdf_interp = griddedInterpolant(agrid, a_cdf, 'linear');

	amin = agrid(1);
	spacing = (amax - amin) / nbins;
	bins = 1:nbins;
	bins = amin + bins * spacing;

	vals = zeros(nbins, 1);
	for ibin = 1:nbins
		bin_start = bins(ibin) - spacing;
		bin_end = bins(ibin);

		P_lt_start = cdf_interp(bin_start);

		if ibin < nbins
			P_lt_end = cdf_interp(bin_end);
		else
			P_lt_end = 1;
		end

		vals(ibin) = P_lt_end - P_lt_start;
	end

	bins = [amin, bins];
end