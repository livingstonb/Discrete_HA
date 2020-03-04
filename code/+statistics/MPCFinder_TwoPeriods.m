classdef MPCFinder_TwoPeriods < handle
	properties
		p;
		grids;
		income;
		model;

		mean_mpc;
	end

	methods
		function obj = MPCFinder_TwoPeriods(p, grids, income, model)
			obj.p = p;
			obj.grids = grids;
			obj.income = income;
			obj.model = model;
		end

		function compute_mpcs(obj)
			xgrid = obj.grids.a.matrix + obj.income.netymat_broadcast;

			sav_baseline = zeros(obj.p.nx_DST,obj.p.nyP,obj.p.nyF,obj.p.nb,obj.p.nyT);
		    for ib = 1:obj.p.nb
		    for iyF = 1:obj.p.nyF
		    for iyP = 1:obj.p.nyP
		        x_iyP_iyF_iyT = xgrid(:,iyP,iyF,ib,:);
		        sav_iyP_iyF_iyT = obj.model.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
		        sav_baseline(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT, [obj.p.nx_DST 1 1 1 obj.p.nyT]);
		    end
		    end
		    end
		    con_baseline = xgrid - sav_baseline;




		    xmpc = obj.grids.a.matrix + obj.income.netymat_broadcast + obj.p.shocks(5);

			sav = zeros(obj.p.nx_DST,obj.p.nyP,obj.p.nyF,obj.p.nb,obj.p.nyT);
		    for ib = 1:obj.p.nb
		    for iyF = 1:obj.p.nyF
		    for iyP = 1:obj.p.nyP
		        x_iyP_iyF_iyT = xmpc(:,iyP,iyF,ib,:);
		        sav_iyP_iyF_iyT = obj.model.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
		        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT, [obj.p.nx_DST 1 1 1 obj.p.nyT]);
		    end
		    end
		    end
		    con_mpc = xmpc - sav;

		    mpcs = (con_mpc - con_baseline) / obj.p.shocks(5);
		    mpcs = reshape(mpcs, [], obj.p.nyT) * obj.income.yTdist;
		    obj.mean_mpc = dot(mpcs(:), obj.model.adist(:));
		end
	end
end