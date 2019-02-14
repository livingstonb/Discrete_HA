function [AYdiff,model] = solve_EGP(beta,p,xgrid,sgrid,agrid_short,...
                                    prefs,income,Iterating,nextmpcshock,prevmodel)
    % This function performs the method of endogenous grid points to find
    % saving and consumption policy functions. It also calls 
    % find_stationary() to compute the stationary distribution over states 
    % via direct methods (rather than simulations) and stores the results 
    % in the 'model' structure.

    agrid = repmat(agrid_short,p.nyP*p.nyF*p.nb,1);
    
    %% CONSTRUCT EXPECTATIONS MATRIX                                     
    betagrid = beta + prefs.betagrid0;
    
    if p.IterateBeta == 1 && p.Display == 1
        msg = sprintf(' %3.3f',betagrid);
        disp([' Trying betagrid =' msg])
    end

    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    if p.r < 0.001
        % Add income so consumption guess is not all zeros
        %extracon = repmat(kron(min(income.netymat,[],2),ones(p.nx,1)),p.nb,1);
        extracon = 0.02 * repmat(xgrid.full(:),p.nb,1);
    else
        extracon = 0;
    end
    con = p.r * repmat(xgrid.full(:),p.nb,1) + extracon;

    % discount factor matrix, 
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if numel(p.risk_aver) > 1
        % IES heterogeneity - nb is number of IES values
        % betagrid is just beta
        betastacked = speye(p.nyP*p.nyF*p.nx*p.nb) * betagrid;
    else
        % beta heterogeneity
        betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));
        betastacked = sparse(diag(betastacked));
    end

    % Expectations operator (conditional on yT)
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if numel(p.risk_aver) == 1
        Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));
        risk_aver_col = kron(p.risk_aver,ones(p.nx*p.nyP*p.nyF,1));
    else
        Emat = kron(prefs.IEStrans,kron(income.ytrans,speye(p.nx)));
    end

    %% EGP ITERATION
    iter = 1;
    cdiff = 1;
    while iter<p.max_iter && cdiff>p.tol_iter
        if iter==1
            conlast = con;
        else
            conlast = conupdate;
        end
        iter = iter + 1;

        % interpolate to get c(x') using c(x)
        
        % c(x)
        conlast = reshape(conlast,[p.nx p.nyP p.nyF p.nb]);
        % c(x')
        c_xp = zeros(p.nx,p.nyP,p.nyF,p.nb,p.nyT);
        
        % x'(s)
        temp_sav = repmat(sgrid.full(:),p.nb,p.nyT);
        temp_inc = repmat(kron(income.netymat,ones(p.nx,1)),p.nb,1);
        xp_s = (1+p.r)*temp_sav + temp_inc + nextmpcshock;
        xp_s = reshape(xp_s,[p.nx p.nyP p.nyF p.nb p.nyT]);

        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            if isempty(prevmodel)
                % usual method of EGP
                coninterp = griddedInterpolant(xgrid.full(:,iyP,iyF),conlast(:,iyP,iyF,ib),'linear');
                xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);
                c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
            else
                % need to compute IMPC(s,t) for s > 1, where IMPC(s,t) is MPC in period t out of period
                % s shock that was learned about in period 1 < s
                % coninterp = griddedInterpolant(xgrid.full(:,iyP,iyF),conlast(:,iyP,iyF,ib),'linear');
                xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);
                %c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
                c_xp(:,iyP,iyF,ib,:) = reshape(prevmodel.coninterp{iyP,iyF,ib}(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
            end
        end
        end
        end
        
        % reshape to take expecation over yT first
        c_xp        = reshape(c_xp,[],p.nyT);
        xp_s	    = reshape(xp_s,[],p.nyT);

        % matrix of next period muc, muc(x',yP',yF)
        mucnext     = prefs.u1(c_xp) - p.temptation/(1+p.temptation) * prefs.u1(xp_s);
        
        % now find muc this period as a function of s:
        % variables defined for each (x,yP,yF,beta) in state space,
        % column vecs of length p.nx*p.nyP*p.nyF*p.nb
        savtaxrate  = (1+p.savtax.*(repmat(sgrid.full(:),p.nb,1)>=p.savtaxthresh));
        mu_consumption = (1+p.r)*betastacked*Emat*(mucnext*income.yTdist);
        mu_bequest     = prefs.beq1(repmat(sgrid.full(:),p.nb,1));
        
        % muc(s(x,yP,yF,beta))
        muc_s = (1-p.dieprob) * mu_consumption ./ savtaxrate...
                                                + p.dieprob * mu_bequest;
                
        % c(s)
        con_s = prefs.u1inv(muc_s);
        
        % x(s) = s + stax + c(s)
        x_s = repmat(sgrid.full(:),p.nb,1)...
                        + p.savtax * max(repmat(sgrid.full(:),p.nb,1)-p.savtaxthresh,0)...
                        + con_s;
        x_s = reshape(x_s,[p.nx p.nyP p.nyF p.nb]);

        % interpolate from x(s) to get s(x)
        sav = zeros(p.nx,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            savinterp = griddedInterpolant(x_s(:,iyP,iyF,ib),sgrid.full(:,iyP,iyF),'linear');
            sav(:,iyP,iyF,ib) = savinterp(xgrid.full(:,iyP,iyF)); 
        end
        end
        end

        % deal with borrowing limit
        sav(sav<p.borrow_lim) = p.borrow_lim;

        % updated consumption function, column vec length of
        % length p.nx*p.nyP*p.nyF*p.nb
        conupdate = repmat(xgrid.full(:),p.nb,1) - sav(:)...
                            - p.savtax * max(sav(:)-p.savtaxthresh,0);

        cdiff = max(abs(conupdate(:)-conlast(:)));
        if mod(iter,50) ==0 && p.Display == 1
            disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end

    end
    
    if cdiff>p.tol_iter
        % EGP did not converge, don't find stationary distribution
        AYdiff = 100000;
        model.EGP_cdiff = cdiff;
        return
    end

    model.sav = sav;
    model.con = reshape(conupdate,[p.nx p.nyP p.nyF p.nb]);
    model.EGP_cdiff = cdiff;
    
    % create interpolants from optimal policy functions
    % and find saving values associated with xvals
    model.savinterp = cell(p.nyP,p.nyF,p.nb);
    model.coninterp = cell(p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        model.savinterp{iyP,iyF,ib} = ...
            griddedInterpolant(xgrid.full(:,iyP,iyF),model.sav(:,iyP,iyF,ib),'linear');
        model.coninterp{iyP,iyF,ib} = ...
            griddedInterpolant(xgrid.full(:,iyP,iyF),model.con(:,iyP,iyF,ib),'linear');    
    end
    end
    end

    %% DISTRIBUTION
    
    
    if Iterating == 1
        % only get distribution over assets
        model.adist = find_stationary_adist(p,model,income,prefs,agrid_short);
    else
        
        [model.adist,model.xdist,model.xvals,model.y_x,model.nety_x,model.statetrans,model.adiff]...
                    = find_stationary_adist(p,model,income,prefs,agrid_short);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP 
            model.sav_x(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(model.xvals(:,iyP,iyF,ib));
        end
        end
        end
        model.sav_x = max(model.sav_x,p.borrow_lim);
        
        % Collapse the asset distribution from (a,yP_lag,yF_lag,beta_lag) to (a,beta_lag) for norisk
        % model, and from (x,yP,yF,beta) to (x,beta)
        if p.nyP>1 && p.nyF>1
            % a
            model.adist_noincrisk =  sum(sum(model.adist,3),2);
            % x
            model.xdist_noincrisk    = sum(sum(model.xdist,3),2);
        elseif (p.nyP>1 && p.nyF==1) || (p.nyP==1 && p.nyF>1)
            model.adist_noincrisk =  sum(model.adist,2);
            model.xdist_noincrisk    = sum(model.xdist,2);
        elseif p.nyP==1 && p.nyF==1
            model.adist_noincrisk = model.adist;
            model.xdist_noincrisk    = model.xdist;
        end
    
        % Policy functions associated with xdist
        model.con_x= model.xvals - model.sav_x - p.savtax*max(model.sav_x-p.savtaxthresh,0);
    end
    
    % mean saving, mean assets
    model.mean_a = model.adist(:)' * agrid(:);
    
    if p.GRIDTEST == 2
        % use simulation results in objective function
        sim = simulate(p,income,model,xgrid,prefs);
        mean_assets = sim.mean_a;
    else
        % use distribution results
        mean_assets = model.mean_a;
    end
           
    if p.Display == 1
        fprintf(' A/Y = %2.3f\n',mean_assets/(income.meany1*p.freq));
    end
    AYdiff = mean_assets/(income.meany1*p.freq) -  p.targetAY;

end