classdef Grid < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties (SetAccess = private)
		x; % cash-on-hand grids
		s; % savings grids
		a; % asset grid
		gtype;
		nx;
		nx_neg;
		nyP;
		nyF;

		i0;
	end

	methods
		function obj = Grid(params, income, gtype)
			obj.gtype = gtype;

			if strcmp(gtype,'EGP')
                obj.nx = params.nx;
                obj.nx_neg = params.nx_neg;
			elseif strcmp(gtype,'DST')
				obj.nx = params.nx_DST;
				obj.nx_neg = params.nx_neg_DST;
			end

			if (params.borrow_lim<0) && (obj.nx_neg==0)
				error("Must have nx_neg > 0 since borrow_lim < 0")
			end

			obj.create_sgrid(params);
			obj.create_xgrid(params, income);
			obj.create_norisk_xgrid(params, income);
            obj.create_agrid(params);

            obj.i0 = obj.nx_neg + 1;
		end

		function obj = create_sgrid(obj, params)
		    soft_constraint = 0;
			if obj.nx_neg > 0
				neg_midpt = (params.borrow_lim + soft_constraint) / 2;
				pts_bottom = round(obj.nx_neg/2);
				pts_top = obj.nx_neg - pts_bottom + 2;

				% Portion of savings grid closest to borrowing limit
				savgrid_neg_low = create_curved_grid(...
					params.borrow_lim, neg_midpt, pts_bottom,...
					params.xgrid_par_neg, false);

				% Portion of savings grid closest to soft constraint
				savgrid_neg_high = create_curved_grid(...
					neg_midpt, soft_constraint, pts_top,...
					params.xgrid_par_neg, true);

				% Combine while deleting repetitions
				savgrid_neg = [savgrid_neg_low; savgrid_neg_high(2:end-1)];
			else
				savgrid_neg = [];
			end

			nx_pos = obj.nx - obj.nx_neg;
			savgrid_pos = create_curved_grid(...
				soft_constraint, params.xmax, nx_pos,...
				params.xgrid_par, false);
			savgrid_pos = obj.enforce_min_spacing(params, savgrid_pos);

			savgrid = [savgrid_neg; savgrid_pos];

		    obj.s.vec = savgrid;
		    obj.s.matrix = repmat(savgrid, [1 params.nyP params.nyF]);
		end

		function obj = create_xgrid(obj, params, income)
			minyT = reshape(min(income.netymat, [], 2), [1 params.nyP params.nyF]);

            R_broadcast = reshape(params.R, [1 1 1 numel(params.R)]);
			xgrid = R_broadcast .* obj.s.matrix + minyT;
			if numel(params.r) == 1
				xgrid = repmat(xgrid, [1 1 1 params.nb]);
			end

			obj.x.matrix = xgrid;
		end

		function obj = create_norisk_xgrid(obj, params, income)
            R_broadcast = reshape(params.R, [1 1 1 numel(params.R)]);
			xmatrix_norisk = R_broadcast .* obj.s.vec + income.meannety1;
			xmatrix_norisk = reshape(xmatrix_norisk, obj.nx, []);

			if numel(params.r) == 1
				xmatrix_norisk = repmat(xmatrix_norisk, [1 params.nb]);
			end

			obj.x.matrix_norisk = xmatrix_norisk;
		end

		function obj = create_agrid(obj, params)
% 			agrid = create_curved_grid(...
% 				params.borrow_lim, params.xmax, obj.nx,...
% 				params.xgrid_par, false);
% 		    obj.a.vec = obj.enforce_min_spacing(params, agrid);
% 		    obj.a.matrix = repmat(obj.a.vec, [1,params.nyP,params.nyF,params.nb]);
		 	obj.a.vec = min(params.R) * obj.s.vec;
		 	obj.a.matrix = min(params.R) * repmat(obj.s.matrix, [1 1 1 params.nb]);
		end

		function grid_adj = enforce_min_spacing(obj, params, gridvec)
			grid_adj = gridvec;
			for ii = 1:numel(gridvec)-1
				if grid_adj(ii+1) - grid_adj(ii) < params.gridspace_min
					grid_adj(ii+1) = grid_adj(ii) + params.gridspace_min;
				else
					break
				end
			end
		end
	end

end

function cgrid = create_curved_grid(lowval, highval, npts,...
	curvature, reversed)
	cgrid = linspace(0, 1, npts)';
    cgrid = cgrid .^ (1/curvature);
    if reversed
    	cgrid = 1 - flip(cgrid);
    end
    cgrid = lowval + (highval - lowval) * cgrid;
end