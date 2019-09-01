classdef EGP_EZ_Solver < handle
    % This class finds the policy functions using the
    % method of endogenous grid points for the case of
    % Epstein-Zin utility
    
	properties (SetAccess = private)
        % parameters and grid
        p;
        grids;
        
		betagrid;
		betastacked;

        % policy and value functions
		con;
		conupdate;
		V;
		Vupdate;
        sav;

		Emat;
		invies_col;
		invies_col_yT;
		risk_aver_col;
		risk_aver_col_yT;
        
        EGP_cdiff = 1e5;
        
        % policy functions of next period cash-on-hand
		c_xp;
		V_xp;

		Vinterp;

		ezval;
	end

	methods
		function obj = EGP_EZ_Solver(betta,p,grids,heterogeneity,income)
            obj.p = p;
            obj.grids = grids;
			obj.betagrid = betta + heterogeneity.betagrid0;

            if obj.p.IterateBeta == 1
		        msg = sprintf(' %3.3f',obj.betagrid);
		        disp([' Trying betagrid =' msg])
            end

		    % initial guess for consumption function, stacked state combinations
		    obj.con = (obj.p.r + 0.002*(obj.p.r<0.001)) * repmat(obj.grids.x.matrix(:),obj.p.nb,1);

		    % initial guess for value function
		    obj.V = obj.con;

		    obj.create_discount_factor_array();
		    obj.create_other_objects(heterogeneity,income);
		end

		function create_discount_factor_array(obj)
			% discount factor matrix, 
		    if (numel(obj.p.invies) > 1) || (numel(obj.p.risk_aver) > 1)
		        obj.betastacked = speye(obj.p.nyP*obj.p.nyF*obj.p.nx*obj.p.nb) * obj.betagrid;
		    else
		        obj.betastacked = kron(obj.betagrid,ones(obj.p.nyP*obj.p.nyF*obj.p.nx,1));
		        obj.betastacked = sparse(diag(obj.betastacked));
		    end
		end

		function create_other_objects(obj,heterogeneity,income)
			% construct xpectations operator (conditional on yT)
		    % construct arrays of invies and riskaver
		    if numel(obj.p.invies) > 1
		        obj.Emat = kron(heterogeneity.ztrans,kron(income.ytrans,speye(obj.p.nx)));
		        obj.invies_col = kron(obj.p.invies',ones(obj.p.nx*obj.p.nyP*obj.p.nyF,1));
		        obj.risk_aver_col = obj.p.risk_aver;
		        obj.invies_col_yT = repmat(obj.invies_col,1,obj.p.nyT);
		    elseif numel(obj.p.risk_aver) > 1
		        obj.Emat = kron(heterogeneity.ztrans,kron(income.ytrans,speye(obj.p.nx)));
		        obj.risk_aver_col = kron(obj.p.risk_aver',ones(obj.p.nx*obj.p.nyP*obj.p.nyF,1));
		        obj.invies_col = obj.p.invies;
		        obj.risk_aver_col_yT = repmat(obj.risk_aver_col,1,obj.p.nyT);
		    else
		        obj.Emat = kron(heterogeneity.betatrans,kron(income.ytrans,speye(obj.p.nx)));
		        obj.risk_aver_col = obj.p.risk_aver;
		        obj.invies_col = obj.p.invies;
		    end
		end

		function solve(obj,income)
			iter = 1;
			while (iter < obj.p.max_iter) && (obj.EGP_cdiff > obj.p.tol_iter)
				obj.iterate_once(income);

				obj.EGP_cdiff = max(abs(obj.conupdate(:)-obj.con(:)));
		        if mod(iter,50) ==0
		            disp([' EGP Iteration ' int2str(iter)...
		            		' max con fn diff is ' num2str(obj.EGP_cdiff)]);
		        end

		        obj.con = obj.conupdate;
		        obj.V = obj.Vupdate;

				iter = iter + 1;
			end
		end

		function iterate_once(obj,income)
			obj.con = reshape(obj.con,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb]);
			obj.V = reshape(obj.V,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb]);

			% store c(x') and V(x')
			obj.update_fns_of_xp(income);

			% matrix of next period muc, muc(x',yP',yF)
			emuc = obj.get_expected_muc(income);

			% current muc(s)
			muc_s = obj.get_current_muc(income,emuc);

			% c(s) by inverting marginal utility
			con_s = muc_s .^ (-1./obj.invies_col);

			% x(s)
			x_s = con_s + repmat(obj.grids.s.matrix(:),obj.p.nb,1)...
                        + obj.p.savtax * max(repmat(obj.grids.s.matrix(:),obj.p.nb,1)-obj.p.savtaxthresh,0);
	        x_s = reshape(x_s,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb]);

	        % s(x) saving policy function
	        obj.sav = obj.get_sav_x_by_interpolating_x_s(x_s);

	        % x'
	        index_to_extend = 1*(obj.p.nyF==1) + 2*(obj.p.nyF>1);
	        xp = obj.p.R * repmat(obj.sav(:),1,obj.p.nyT) ... 
	                    + repmat(kron(income.netymat,ones(obj.p.nx,1)),obj.p.nb,index_to_extend);
	        xp = reshape(xp,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb obj.p.nyT]);

	        % update consumption
	        obj.conupdate = repmat(obj.grids.x.matrix,[1 1 1 obj.p.nb]) ...
	        	- obj.sav - obj.p.savtax * max(obj.sav-obj.p.savtaxthresh,0);

	        % compute E[V(x)^(1-riskaver)]^(1/(1-riskaver))
	        obj.update_ezval(income,xp);

	        obj.update_value_fn();
		end

		function update_fns_of_xp(obj,income)
			obj.c_xp = zeros(obj.p.nx,obj.p.nyP,obj.p.nyF,obj.p.nb,obj.p.nyT);
			obj.V_xp = zeros(obj.p.nx,obj.p.nyP,obj.p.nyF,obj.p.nb,obj.p.nyT);

			temp_sav = repmat(obj.grids.s.matrix(:),obj.p.nb,obj.p.nyT);
	        temp_inc = repmat(kron(income.netymat,ones(obj.p.nx,1)),obj.p.nb,1);
	        xp_s = (1+obj.p.r)*temp_sav + temp_inc;
	        xp_s = reshape(xp_s,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb obj.p.nyT]);

            for ib  = 1:obj.p.nb
            for iyF = 1:obj.p.nyF
            for iyP = 1:obj.p.nyP
	            xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);
	            coninterp = griddedInterpolant(obj.grids.x.matrix(:,iyP,iyF),obj.con(:,iyP,iyF,ib),'linear');
	            obj.c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,obj.p.nyT);
	            obj.Vinterp{iyP,iyF,ib} = griddedInterpolant(obj.grids.x.matrix(:,iyP,iyF),obj.V(:,iyP,iyF,ib),'linear');
	            obj.V_xp(:,iyP,iyF,ib,:) = reshape(obj.Vinterp{iyP,iyF,ib}(xp_s_ib_iyF_iyP(:)),[],1,1,1,obj.p.nyT);
            end
            end
            end

	        obj.c_xp = reshape(obj.c_xp,[],obj.p.nyT);
	        obj.V_xp = reshape(obj.V_xp,[],obj.p.nyT);
		end

		function emuc = get_expected_muc(obj,income)
			% nexts period's muc(x',yP',yF)
	        if numel(obj.p.invies) > 1
	            mucnext = obj.c_xp.^(-obj.invies_col_yT) ...
	            	.* obj.V_xp.^(obj.invies_col_yT-obj.p.risk_aver);
	        elseif numel(obj.p.risk_aver) > 1
	            mucnext = obj.c_xp.^(-obj.p.invies) .* obj.V_xp.^(obj.p.invies-obj.risk_aver_col_yT);
	        else
	            mucnext = obj.c_xp.^(-obj.p.invies) .* obj.V_xp.^(obj.p.invies-obj.p.risk_aver);
	        end
	        
	        % expected muc
	        savtaxrate  = (1+obj.p.savtax.*(repmat(obj.grids.s.matrix(:),obj.p.nb,1)>=obj.p.savtaxthresh));
	        mu_cons = (1+obj.p.r)*obj.betastacked*obj.Emat*mucnext*income.yTdist ./ savtaxrate;
	        mu_bequest = aux.utility_bequests1(obj.p.bequest_curv,obj.p.bequest_weight,...
	    		obj.p.bequest_luxury,repmat(obj.grids.s.matrix(:),obj.p.nb,1));
	        emuc = (1-obj.p.dieprob) * mu_cons + obj.p.dieprob * mu_bequest;
		end

		function muc_s = get_current_muc(obj,income,emuc)
            if numel(obj.p.risk_aver) > 1
	            ezvalnext_ra_equal1 = exp(obj.Emat * log(obj.V_xp) * income.yTdist)...
	                .^ (obj.risk_aver_col-obj.invies_col);

	            ezvalnext_ra_nequal1 = (obj.Emat * obj.V_xp.^(1-obj.risk_aver_col)...
	            	* income.yTdist)...
	            	.^ ((obj.risk_aver_col-obj.invies_col)./(1-obj.risk_aver_col));

	            ezvalnext = zeros(obj.p.nx*obj.p.nyP*obj.p.nyF*obj.p.nb,1);
	            ezvalnext(obj.risk_aver_col==1) = ezvalnext_ra_equal1(obj.risk_aver_col==1);
	            ezvalnext(obj.risk_aver_col~=1) = ezvalnext_ra_nequal1(obj.risk_aver_col~=1);
            else
                if obj.p.risk_aver == 1
	                ezvalnext = exp(obj.Emat * log(obj.V_xp) * income.yTdist)...
	                    .^ (obj.risk_aver_col - obj.invies_col);
	            else
	                ezvalnext = (obj.Emat * obj.V_xp .^ (1-obj.risk_aver_col) * income.yTdist)...
	                    .^ ((obj.risk_aver_col - obj.invies_col)./(1 - obj.risk_aver_col));
                end
            end
	        
	        muc_s = emuc .* ezvalnext;
	    end

	    function sav = get_sav_x_by_interpolating_x_s(obj,x_s)
	    	% interpolate from x(s) to get s(x)
	        sav = zeros(obj.p.nx,obj.p.nyP,obj.p.nyF,obj.p.nb);
	        for ib  = 1:obj.p.nb
	        for iyF = 1:obj.p.nyF
	        for iyP = 1:obj.p.nyP
	            savinterp = griddedInterpolant(x_s(:,iyP,iyF,ib),obj.grids.s.matrix(:,iyP,iyF),'linear');
	            sav(:,iyP,iyF,ib) = savinterp(obj.grids.x.matrix(:,iyP,iyF)); 
	        end
	        end
	        end
	        sav = max(sav,obj.p.borrow_lim);
	    end

	    function update_ezval(obj,income,xp)
	    	% interpolate adjusted expected value function on x grid
	        ezval_integrand = zeros(obj.p.nx,obj.p.nyP,obj.p.nyF,obj.p.nb,obj.p.nyT);
	        for ib = 1:numel(obj.p.risk_aver)
	        for iyF = 1:obj.p.nyF
	        for iyP = 1:obj.p.nyP
	            xp_iyP_iyF_ib = xp(:,iyP,iyF,ib,:);
	            temp_iyP_iyF_ib = obj.Vinterp{iyP,iyF,ib}(xp_iyP_iyF_ib(:)) .^ (1-obj.p.risk_aver(ib));
	            ezval_integrand(:,iyP,iyF,ib,:) = reshape(temp_iyP_iyF_ib,[obj.p.nx 1 1 1 obj.p.nyT]);
	        end
	        end
	        end

	        % Take expectation over yT
	        ezval_integrand = reshape(ezval_integrand,[],obj.p.nyT) * income.yTdist;
	        % Take expectation over (yP,yF,beta)
	        obj.ezval = obj.Emat * ezval_integrand;

	        if numel(obj.p.risk_aver) > 1
	            ezval_ra_equal1 = (obj.risk_aver_col==1) .* exp(obj.ezval);
	            ezval_ra_nequal1 = (obj.risk_aver_col~=1) ...
	            	.* obj.ezval .^ (1./(1-obj.risk_aver_col));

	            obj.ezval = zeros(obj.p.nx*obj.p.nyP*obj.p.nyF*obj.p.nb,1);
	            obj.ezval(obj.risk_aver_col==1) = ezval_ra_equal1(obj.risk_aver_col==1);
	            obj.ezval(obj.risk_aver_col~=1) = ezval_ra_nequal1(obj.risk_aver_col~=1);
	        else
	            if obj.p.risk_aver == 1
	                obj.ezval = exp(obj.ezval);
	            else
	                obj.ezval = obj.ezval .^(1./(1-obj.p.risk_aver));
	            end
	        end

	        obj.ezval = reshape(obj.ezval,obj.p.nx,obj.p.nyP,obj.p.nyF,obj.p.nb);
	    end

	    function update_value_fn(obj)
	    	obj.Vupdate = zeros(obj.p.nx,obj.p.nyP,obj.p.nyF,obj.p.nb);
	        for ib = 1:obj.p.nb
	            if numel(obj.p.invies) == 1

	                if numel(obj.p.risk_aver) == 1
	                    ibeta = ib; % possible beta heterogeneity
	                else
	                    ibeta = 1;
	                end

	                if obj.p.invies == 1
	                    obj.Vupdate(:,:,:,ib) = obj.conupdate(:,:,:,ib) .^ (1-obj.betagrid(ibeta)) ...
	                    	.* obj.ezval(:,:,:,ib) .^ obj.betagrid(ibeta);
	                else
	                	obj.Vupdate(:,:,:,ib) = (1-obj.betagrid(ibeta)) * obj.conupdate(:,:,:,ib) .^ (1-obj.p.invies) ...
	                                    + obj.betagrid(ibeta) * obj.ezval(:,:,:,ib) .^ (1-obj.p.invies);
	                    obj.Vupdate(:,:,:,ib) = obj.Vupdate(:,:,:,ib) .^ (1/(1-obj.p.invies));
	                end
	            elseif numel(obj.p.invies) > 1
	                if obj.p.invies(ib) == 1
	                    obj.Vupdate(:,:,:,ib) = obj.conupdate(:,:,:,ib) .^ (1-obj.betagrid) ...
	                    	.* obj.ezval(:,:,:,ib) .^ obj.betagrid;
	                else
	                    obj.Vupdate(:,:,:,ib) = (1-obj.betagrid) * obj.conupdate(:,:,:,ib) .^ (1-obj.p.invies(ib)) ...
	                                    + obj.betagrid * obj.ezval(:,:,:,ib) .^ (1-obj.p.invies(ib));
	                    obj.Vupdate(:,:,:,ib) = obj.Vupdate(:,:,:,ib) .^ (1/(1-obj.p.invies(ib)));
	                end
	            end
	        end

	        assert(sum(~isfinite(obj.Vupdate(:)))==0)
	        assert(all(obj.Vupdate(:)>=0))
	    end

	    function model = return_model(obj)
	    	model = struct();
	    	model.sav = obj.sav;
		    model.con = reshape(obj.con,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb]);
		    model.V = reshape(obj.V,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb]);
		    
		    % create interpolants from optimal policy functions
		    % and find saving values associated with xvals
		    model.savinterp = cell(obj.p.nyP,obj.p.nyF,obj.p.nb);
		    model.coninterp = cell(obj.p.nyP,obj.p.nyF,obj.p.nb);
		    for ib = 1:obj.p.nb
		    for iyF = 1:obj.p.nyF
		    for iyP = 1:obj.p.nyP
		        model.savinterp{iyP,iyF,ib} = ...
		            griddedInterpolant(obj.grids.x.matrix(:,iyP,iyF),model.sav(:,iyP,iyF,ib),'linear');
		        model.coninterp{iyP,iyF,ib} = ...
		            griddedInterpolant(obj.grids.x.matrix(:,iyP,iyF),model.con(:,iyP,iyF,ib),'linear');    
		    end
		    end
            end
            
            model.EGP_cdiff = obj.EGP_cdiff;
	    end
	end

end