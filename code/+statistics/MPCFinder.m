classdef MPCFinder < handle
	% This class is used for computing MPCs after solving for
	% the policy functions and stationary distribution.

	% Results are storied in 'mpcs', which is filled with NaNs
	% after instantiation. Next, the solve() method is called
	% and the results can be retrieved.

	% 'basemodel' is a structure containing the policy functions
	% and the stationary distribution of the model

	% 'models' is a cell array indexed by (1) the shock index,
	% (2) the shock period, and (3) the current period. This
	% array contains the policy functions based on expectation
	% of a future shock.

	properties (SetAccess = private)
		Nstates; % total number of states
		r_mat; % matrix of interest rates, if more than 1

		income;

		con_baseline; % baseline consumption
		con_baseline_yT; % baseline consumption before expectation wrt yT
		basemodel; % baseline model (policy fns, etc)
		models; %  models for other shock timings

		xgrid_yT; % cash-on-hand, includes a yT dimension

		fspace; % object used to interpolate back onto agrid

		mpcs = struct(); % results
		loan;
		loss_in_2_years;
	end

	methods
		function obj = MPCFinder(p,income,grids,basemodel,models)
			obj.Nstates = p.nx_DST*p.nyP*p.nyF*p.nb;
			obj.basemodel = basemodel;
			obj.models = models;
			obj.fspace = fundef({'spli',grids.a.vec,0,1});
			obj.income = income;

			obj.xgrid_yT = grids.a.vec + income.netymat_broadcast;

			if numel(p.r) > 1
		        obj.r_mat = reshape(r_col,[1 1 1 numel(p.r)]);
		    else
		        obj.r_mat = p.r;
		    end

		    for ishock = 1:6
		    	% mean mpc in period t out of period s shock
		    	obj.mpcs(ishock).avg_s_t = NaN(5,5);
		    	% mpcs over states for shock in period 1
		    	obj.mpcs(ishock).mpcs_1_t = cell(1,4);
                % fraction of responders
                obj.mpcs(ishock).mpc_pos = NaN(5,5);
                obj.mpcs(ishock).mpc_neg = NaN(5,5);
                obj.mpcs(ishock).mpc0 = NaN(5,5);
                % conditional mean
                obj.mpcs(ishock).avg_s_t_condl = NaN(5,5);

		    	obj.loan = NaN;
		    	obj.loss_in_2_years = NaN;
	    	end
		end

		function solve(obj,p,grids)
			% This function calls the appropriate methods to find
			% the MPCs.

			obj.get_baseline_consumption(p);

			for ishock = 1:6
				if (p.mpcshocks_after_period1 == 0) || (p.EpsteinZin == 1)
					shockperiods = 1;
				else
					shockperiods = [1 2 5];
				end

				for shockperiod = shockperiods
					loan = 0;
					obj.computeMPCs(p,grids,ishock,shockperiod,loan);
				end
            end

            if p.mpcshocks_after_period1 > 0
                % $500 loss in 2 years
                obj.computeMPCs(p,grids,1,9,0);

                % $5000 loan for one year
                obj.computeMPCs(p,grids,3,5,0.081);
            end

			obj.compute_cumulative_mpcs();
		end

		function get_baseline_consumption(obj,p)
			% Each (a,yP,yF) is associated with nyT possible x values, create this
		    % grid here
		    obj.con_baseline_yT = obj.get_policy(p,obj.xgrid_yT,obj.basemodel);
		    % take expectation wrt yT
		    obj.con_baseline = reshape(obj.con_baseline_yT,[],p.nyT) * obj.income.yTdist;
		end

		function con = get_policy(obj,p,x_mpc,model)
			% Computes consumption policy function after taking expectation to get
		    % rid of yT dependence
		    sav = zeros(p.nx_DST,p.nyP,p.nyF,p.nb,p.nyT);
		    for ib = 1:p.nb
		    for iyF = 1:p.nyF
		    for iyP = 1:p.nyP
		        x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,1,:);
		        sav_iyP_iyF_iyT = model.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
		        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,[p.nx_DST 1 1 1 p.nyT]);
		    end
		    end
		    end
		    sav = max(sav,p.borrow_lim);
		    x_mpc = reshape(x_mpc,[p.nx_DST p.nyP p.nyF 1 p.nyT]);
		    con = x_mpc - sav - p.savtax * max(sav-p.savtaxthresh,0);
		end

		function computeMPCs(obj,p,grids,ishock,shockperiod,loan)
			shock = p.shocks(ishock);

			for it = 1:shockperiod
				% transition matrix from period 1 to period t
				if it == 1
					trans_1_t = speye(obj.Nstates);
				else
					previous_period = it - 1;
					one_period_transition = obj.transition_matrix_given_t_s(...
						p,shockperiod,previous_period,ishock);
 					trans_1_t =  one_period_transition * trans_1_t;
				end

				% cash-on-hand
				if it == shockperiod
					x_mpc = obj.xgrid_yT + shock;
				elseif it == 1
					x_mpc = obj.xgrid_yT + loan;
				else
					x_mpc = obj.xgrid_yT;
				end

				if shock < 0
	            	% record which states are pushed below asset grid after negative shock
	                % bring back up to asset grid for interpolation
	                below_xgrid = false(size(x_mpc));
	                for iyT = 1:p.nyT
	                    below_xgrid(:,:,:,:,iyT) = x_mpc(:,:,:,:,iyT) < grids.x.matrix(1,:,:);

	                    x_mpc(:,:,:,:,iyT) = ~below_xgrid(:,:,:,:,iyT) .* x_mpc(:,:,:,:,iyT)...
	                                        + below_xgrid(:,:,:,:,iyT) .* grids.x.matrix(1,:,:);
	                end
	                below_xgrid = reshape(below_xgrid,[p.nx_DST p.nyP p.nyF 1 p.nyT]);
	            end

	            % consumption choice given the shock
	            con = obj.get_policy(p,x_mpc,obj.models{ishock,shockperiod,it});

	            if (shock < 0) && (it == shockperiod)
	                % make consumption for cases pushed below xgrid equal to consumption
	                % at bottom of xgrid - the amount borrowed
	                x_before_shock = reshape(grids.x.matrix,[p.nx_DST p.nyP p.nyF]);
	                x_minus_xmin = x_before_shock - grids.x.matrix(1,:,:);
	            	con = ~below_xgrid .* con ...
	            		+ below_xgrid .* (obj.con_baseline_yT(1,:,:,:,:)...
	            							+ shock + x_minus_xmin);
	            end

	            % expectation over yT
	            con = reshape(con,[],p.nyT) * obj.income.yTdist;

	            % now compute IMPC(s,t)
	            mpcs = ( trans_1_t * con - obj.basemodel.statetrans^(it-1) * obj.con_baseline) / shock;

	            if loan > 0
	            	obj.loan = obj.basemodel.adist(:)' * mpcs(:);
	            elseif shockperiod <= 5
	            	obj.mpcs(ishock).avg_s_t(shockperiod,it) = obj.basemodel.adist(:)' * mpcs(:);
                    
                    loc_pos = mpcs(:) > 0;
                    dist_vec = obj.basemodel.adist(:);
                    if sum(dist_vec(loc_pos)) > 0
                        obj.mpcs(ishock).avg_s_t_condl(shockperiod,it) = dist_vec(loc_pos)' * mpcs(loc_pos) / sum(dist_vec(loc_pos));
                    end
                    obj.mpcs(ishock).mpc_pos = sum(dist_vec(loc_pos));
                    obj.mpcs(ishock).mpc_neg = sum(dist_vec(mpcs(:)<0));
                    obj.mpcs(ishock).mpc0 = sum(dist_vec(mpcs(:)==0));
	            else
	            	obj.loss_in_2_years = obj.basemodel.adist(:)' * mpcs(:);
	            end
	            
	            if (shockperiod == 1) && (it >= 1 && it <= 4)
	            	% store is = 1 mpcs for decompositions
	                obj.mpcs(ishock).mpcs_1_t{it} = mpcs;
	            end
			end

			obj.computeMPCs_periods_after_shock(p,grids,ishock,...
				shockperiod,trans_1_t);
		end

		function computeMPCs_periods_after_shock(obj,p,grids,...
			ishock,shockperiod,trans_1_t)
			shock = p.shocks(ishock);

			% transition probabilities from it = is to it = is + 1
			trans_1_t = trans_1_t * obj.transition_matrix_given_t_s(p,shockperiod,shockperiod,ishock);

	        RHScon = obj.basemodel.statetrans^shockperiod * obj.con_baseline(:);
	        LHScon = obj.con_baseline(:);

	        for it = shockperiod+1:5 % it > shockperiod case, policy fcns stay the same in this region
	            mpcs = (trans_1_t * LHScon - RHScon) / shock;

	            % state transitions are as usual post-shock
	            obj.mpcs(ishock).avg_s_t(shockperiod,it) = obj.basemodel.adist(:)' * mpcs(:);
	            RHScon = obj.basemodel.statetrans * RHScon;
	            LHScon = obj.basemodel.statetrans * LHScon;
	            
	            if (shockperiod == 1) && (it >= 1 && it <= 4)
	                obj.mpcs(ishock).mpcs_1_t{it} = mpcs;
	            end
	        end
		end

		function compute_cumulative_mpcs(obj)
			for ishock = 1:6
				obj.mpcs(ishock).avg_1_1to4 = sum(obj.mpcs(ishock).avg_s_t(1,1:4));
				obj.mpcs(ishock).avg_5_1to4 = sum(obj.mpcs(ishock).avg_s_t(5,1:4));
			end
		end

		function transition = transition_matrix_given_t_s(...
			obj,p,is,ii,ishock)
			% Computes the transition matrix between t=ii and 
			% t=ii + 1 given shock in period 'is'

			% shocked cash-on-hand
            if is == ii
                shock = p.shocks(ishock);
            else
                shock = 0;
            end
			x_mpc = obj.xgrid_yT + shock;

			% get saving policy function
			sav = zeros(p.nx_DST,p.nyP,p.nyF,p.nb,p.nyT);
			reshape_vec = [p.nx_DST 1 1 1 p.nyT];
			for ib = 1:p.nb
			for iyF = 1:p.nyF
			for iyP = 1:p.nyP
				x_iyP_iyF_iyT = reshape(x_mpc(:,iyP,iyF,1,:),[],1);

				if (shock < 0) && (ii == is)
		        	below_xgrid = x_iyP_iyF_iyT < min(obj.xgrid_yT(1,iyP,iyF,1,:));
		            x_iyP_iyF_iyT(below_xgrid) = min(obj.xgrid_yT(1,iyP,iyF,1,:));
                end
		        
                sav_iyP_iyF_iyT = obj.models{ishock,is,ii}.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT);
                
		        if (shock < 0) && (ii == is)
		            sav_iyP_iyF_iyT(below_xgrid) = 0;
		        end

		        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,reshape_vec);
			end
			end
			end

			sav = max(sav,p.borrow_lim);

			% next period's assets conditional on living
			aprime_live = (1+repmat(obj.r_mat,[1 1 1 1 p.nyT])) .* sav;

			% interpolate next period's assets back onto asset grid
			asset_interp = funbas(obj.fspace,aprime_live(:));
		    asset_interp = reshape(asset_interp,obj.Nstates,p.nx_DST*p.nyT)...
		    	* kron(speye(p.nx_DST),obj.income.yTdist);
		    if p.Bequests == 1
		        interp_death = asset_interp;
		    else
		        interp_death = sparse(obj.Nstates,p.nx_DST);
		        interp_death(:,1) = 1;
            end
            
            ytrans_live_long = kron(obj.income.ytrans_live, ones(p.nx_DST,1));
            ytrans_death_long = kron(obj.income.ytrans_death, ones(p.nx_DST,1));

		    % now construct transition matrix
		    transition = sparse(obj.Nstates,obj.Nstates);
		    for col = 1:p.nyP*p.nyF*p.nb
		        newblock_live = ytrans_live_long(:,col) .* asset_interp;
		        newblock_death = ytrans_death_long(:,col) .* interp_death;
		        transition(:,p.nx_DST*(col-1)+1:p.nx_DST*col) = ...
		        	(1-p.dieprob)*newblock_live + p.dieprob*newblock_death;
		    end
		end
	end

end