function [AYdiff,model] = solve_EGP(beta,p,grids,gridsKFE,prefs,income,nextmpcshock,prevmodel)
    % This function performs the method of endogenous grid points to find
    % saving and consumption policy functions. It also calls 
    % find_stationary() to compute the stationary distribution over states 
    % via direct methods (rather than simulations) and stores the results 
    % in the 'model' structure.

    
    %% CONSTRUCT EXPECTATIONS MATRIX                                     
    betagrid = beta + prefs.betagrid0;
    
    if p.IterateBeta == 1 && p.Display == 1
        msg = sprintf(' %3.3f',betagrid);
        disp([' Trying betagrid =' msg])
    end

    % Expectations operator (conditional on yT)
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if numel(p.r) > 1
        Emat = kron(prefs.rtrans,kron(income.ytrans,speye(p.nx)));
        r_col = kron(p.r',ones(p.nx*p.nyP*p.nyF,1));
        r_mat = reshape(r_col,[p.nx,p.nyP,p.nyF,numel(p.r)]);
    elseif numel(p.risk_aver) > 1
        Emat = kron(prefs.IEStrans,kron(income.ytrans,speye(p.nx)));
        risk_aver_col = kron(p.risk_aver',ones(p.nx*p.nyP*p.nyF,1));
        r_mat = p.r;
    else
        Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));
        r_mat = p.r;
    end

    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    if p.temptation > 0.05
        extra = 0.5;
    else
        extra = 0;
    end
    
    con = (r_mat(:) + 0.002 * (r_mat(:)<0.001) + extra) .* repmat(grids.x.matrix(:),p.nb,1);

    % discount factor matrix, 
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if (numel(p.risk_aver) > 1) || (numel(p.r) > 1)
        % IES heterogeneity or returns heterogeneity - nb is number of IES or r values
        % betagrid is just beta
        betastacked = speye(p.nyP*p.nyF*p.nx*p.nb) * betagrid;
    else
        % beta heterogeneity
        betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));
        betastacked = sparse(diag(betastacked));
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
        temp_sav = repmat(grids.s.matrix(:),p.nb,p.nyT);
        temp_sav = reshape(temp_sav,[p.nx p.nyP p.nyF p.nb p.nyT]);
        index_to_extend = 1*(p.nyF==1) + 2*(p.nyF>1);
        repscheme = ones(1,2);
        repscheme(index_to_extend) = p.nb;
        temp_inc = repmat(kron(income.netymat,ones(p.nx,1)),repscheme);
        temp_inc = reshape(temp_inc,[p.nx p.nyP p.nyF p.nb p.nyT]);
        xp_s = (1+r_mat) .* temp_sav + temp_inc + nextmpcshock;

        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            if isempty(prevmodel)
                % usual method of EGP
                coninterp = griddedInterpolant(grids.x.matrix(:,iyP,iyF),conlast(:,iyP,iyF,ib),'linear');
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
        c_xp = reshape(c_xp,[],p.nyT);
        xp_s = reshape(xp_s,[],p.nyT);

        % matrix of next period muc, muc(x',yP',yF)
        if numel(p.risk_aver) > 1
            risk_aver_col_yT = repmat(risk_aver_col,1,p.nyT);
            mucnext = prefs.u1(risk_aver_col_yT,c_xp)...
                - p.temptation/(1+p.temptation) * prefs.u1(risk_aver_col_yT,xp_s);
        else
            mucnext = prefs.u1(c_xp) - p.temptation/(1+p.temptation) * prefs.u1(xp_s);

        end
            
        
        % now find muc this period as a function of s:
        % variables defined for each (x,yP,yF,beta) in state space,
        % column vecs of length p.nx*p.nyP*p.nyF*p.nb
        savtaxrate  = (1+p.savtax.*(repmat(grids.s.matrix(:),p.nb,1)>=p.savtaxthresh));
        mu_consumption = (1+r_mat(:)).*betastacked*Emat*(mucnext*income.yTdist);
        mu_bequest     = prefs.beq1(repmat(grids.s.matrix(:),p.nb,1));
        
        % muc(s(x,yP,yF,beta))
        muc_s = (1-p.dieprob) * mu_consumption ./ savtaxrate...
                                                + p.dieprob * mu_bequest;
                
        % c(s)
        if numel(p.risk_aver) == 1
            con_s = prefs.u1inv(muc_s);
        else
            con_s = prefs.u1inv(risk_aver_col,muc_s);
        end
        
        % x(s) = s + stax + c(s)
        x_s = repmat(grids.s.matrix(:),p.nb,1)...
                        + p.savtax * max(repmat(grids.s.matrix(:),p.nb,1)-p.savtaxthresh,0)...
                        + con_s;
        x_s = reshape(x_s,[p.nx p.nyP p.nyF p.nb]);

        % interpolate from x(s) to get s(x)
        sav = zeros(p.nx,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            savinterp = griddedInterpolant(x_s(:,iyP,iyF,ib),grids.s.matrix(:,iyP,iyF),'linear');
            sav(:,iyP,iyF,ib) = savinterp(grids.x.matrix(:,iyP,iyF)); 
        end
        end
        end

        % deal with borrowing limit
        sav(sav<p.borrow_lim) = p.borrow_lim;

        % updated consumption function, column vec length of
        % length p.nx*p.nyP*p.nyF*p.nb
        conupdate = repmat(grids.x.matrix(:),p.nb,1) - sav(:)...
                            - p.savtax * max(sav(:)-p.savtaxthresh,0);

        cdiff = max(abs(conupdate(:)-conlast(:)));
        if mod(iter,50) ==0 && p.Display == 1
            disp(['  EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
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
            griddedInterpolant(grids.x.matrix(:,iyP,iyF),model.sav(:,iyP,iyF,ib),'linear');
        model.coninterp{iyP,iyF,ib} = ...
            griddedInterpolant(grids.x.matrix(:,iyP,iyF),model.con(:,iyP,iyF,ib),'linear');    
    end
    end
    end

    %% DISTRIBUTION
    
    model = find_stationary_adist(p,model,income,prefs,gridsKFE);
    
    % get saving policy function defined on xgrid
    model.sav_x = zeros(p.nx_KFE*p.nyT,p.nyP,p.nyF,p.nb);
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
    
    % mean saving, mean assets
	model.mean_a = model.adist(:)' * gridsKFE.a.matrix(:);
    
    if p.GRIDTEST == 2
        % use simulation results in objective function
        sim = simulate(p,income,model,gridsKFE.x.matrix,prefs);
        mean_assets = sim.mean_a;
    else
        % use distribution results
        mean_assets = model.mean_a;
    end
           
    if p.Display == 1
        fprintf(' A/Y = %2.5f\n',mean_assets/(income.meany1*p.freq));
    end
    AYdiff = mean_assets/(income.meany1*p.freq) -  p.targetAY;

end