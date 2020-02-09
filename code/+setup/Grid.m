classdef Grid < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties (SetAccess = private)
		x; % cash-on-hand grids
		s; % savings grids
		a; % asset grid
		gtype;
		nx;
		nyP;
		nyF;
	end

	methods
		function obj = Grid(params, income, gtype)
			obj.gtype = gtype;

			if strcmp(gtype,'EGP')
                obj.nx = params.nx;
			elseif strcmp(gtype,'DST')
				obj.nx = params.nx_DST;
			end

			obj.create_sgrid(params);
			obj.create_xgrid(params, income);
			obj.create_norisk_xgrid(params, income);
            obj.create_agrid(params);
		end

		function obj = create_sgrid(obj, params)
			savgrid = create_curved_grid(...
				params.borrow_lim, params.xmax, obj.nx, params.xgrid_par);
		    savgrid = obj.enforce_min_spacing(params, savgrid);

		    obj.s.vec = savgrid;
		    obj.s.matrix = repmat(savgrid, [1 params.nyP params.nyF]);
		end

		function obj = create_xgrid(obj, params, income)
			minyT = kron(min(income.netymat, [], 2), ones(obj.nx,1));
			xgrid = min(params.R) * obj.s.matrix(:) + minyT;

		    obj.x.matrix = reshape(xgrid, [obj.nx params.nyP params.nyF]);
		end

		function obj = create_norisk_xgrid(obj, params, income)
			obj.x.vec_norisk  = obj.s.vec + income.meannety1;
		end

		function obj = create_agrid(obj, params)
			agrid = create_curved_grid(...
				params.borrow_lim, params.xmax, obj.nx, params.xgrid_par);
		    obj.a.vec = obj.enforce_min_spacing(params, agrid);
		    obj.a.matrix = repmat(obj.a.vec, [1,params.nyP,params.nyF,params.nb]);
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

function cgrid = create_curved_grid(lowval, highval, npts, curvature)
	cgrid = linspace(0, 1, npts)';
    cgrid = cgrid .^ (1/curvature);
    cgrid = lowval + (highval - lowval) * cgrid;
end