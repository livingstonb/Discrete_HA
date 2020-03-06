classdef VFI < handle

	properties
		p;
		grids;
		income;

		ra_bc;
		beta_bc;

		c_policy;
		c_interp;

		ss_dims;

		max_iters = 1e3;
		tol = 1e-9;
	end

	methods
		function obj = VFI(p, grids, income, model, heterogeneity)
			obj.p = p;
			obj.grids = grids;
			obj.income = income;

			obj.c_policy = model.con;
			obj.c_interp = model.coninterp;

			obj.R_bc = heterogeneity.R_broadcast;
			obj.ra_bc = heterogeneity.risk_aver_broadcast;
			obj.beta_bc = heterogeneity.betagrid_broadcast;

			obj.ss_dims = [p.nx, p.nyP, p.nyF, p.nb];
		end

		function value_function = solve(obj)
			Emat = kron(obj.income.ytrans_live, speye(obj.p.nx));

			value_function = obj.make_guess();
			ydist_bc = reshape(obj.income.yTdist, [1, 1, 1, 1, obj.p.nyT]);

			util = aux.utility(obj.ra_bc, obj.c_policy);
			sav = obj.grids.x.matrix - obj.c_policy;
			xprime = obj.R_bc .* sav + obj.income.netymat_broadcast;

			dif = 1e5;
			it = 0;
			while (it <= obj.max_iters) && (dif >= tol)
				it = it + 1;

				val_next = zeros(obj.p.nx, obj.p.nyP, obj.p.nyF, obj.p.nb, obj.p.nyT);
				for ib = 1:obj.p.nb
				for iyF = 1:obj.p.nyF
				for iyP = 1:obj.p.nyP
					val_interp = griddedInterpolant(...
						obj.grids.x.matrix(:,iyP,iyF,ib),...
						value_function(:,iyP,iyF,ib), 'linear');

					xp = xprime(:,iyP,iyF,ib,:);
					val_next(:,iyP,iyF,ib,:) = reshape(val_interp(xp(:)),...
						[obj.p.nx, 1, 1, 1, obj.p.nyT]);
				end
				end
				end

				% Take expectation over yT
				val_next = val_next * ydist_bc;
				
				% Take expectation over remaining random variables
				Eval = reshape(Emat * val_next(:), obj.ss_dims);

				value_update = util + obj.beta_bc .* (1 - obj.p.dieprob) .* Eval;
				dif = max(abs(value_function(:) - value_update(:)));
				value_function = value_update;
			end
		end

		function Vguess = make_guess(obj)
			const = (1 - obj.p.dieprob) * min(obj.beta_bc, 0.995);
			const = const ./ (1 - const);
			Vguess = const .* aux.utility(obj.ra_bc, obj.c_policy);
		end
	end

end