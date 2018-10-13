function [AYdiff,model] = solve_EGP(beta,p,xgrid,sgrid,prefs,income)
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

    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    con = p.r * repmat(xgrid.full(:),p.nb,1);

    % discount factor matrix, 
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));
    betastacked = sparse(diag(betastacked));

    % Expectations operator (conditional on yT)
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));

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
        conlast = reshape(conlast,[p.ns p.nyP p.nyF p.nb]);
        % c(x')
        c_xp = zeros(p.ns,p.nyP,p.nyF,p.nb,p.nyT);
        
        % x'(s)
        temp_sav = repmat(sgrid.full(:),p.nb,p.nyT);
        temp_inc = repmat(kron(income.netymat,ones(p.ns,1)),p.nb,1);
        xp_s = (1+p.r)*temp_sav + temp_inc;
        xp_s = reshape(xp_s,[p.ns p.nyP p.nyF p.nb p.nyT]);

        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            coninterp = griddedInterpolant(xgrid.full(:,iyP,iyF),conlast(:,iyP,iyF,ib),'linear');
            xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);
            c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
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
        muc_savtaxrate  = (1+p.savtax.*(repmat(sgrid.full(:),p.nb,1)>=p.savtaxthresh));
        muc_consumption = (1+p.r)*betastacked*Emat*(mucnext*income.yTdist);
        muc_bequest     = prefs.beq1(repmat(sgrid.full(:),p.nb,1));
        
        % muc(s(x,yP,yF,beta))
        muc_s           = (1-p.dieprob) * muc_consumption ./ muc_savtaxrate...
                                                + p.dieprob * muc_bequest;
                
        % c(s)
        con_s       = prefs.u1inv(muc_s);
        
        % x(s) = s + stax + c(s)
        x_s         = repmat(sgrid.full(:),p.nb,1)...
                        + p.savtax * max(repmat(sgrid.full(:),p.nb,1)-p.savtaxthresh,0)...
                        + con_s;
        x_s         = reshape(x_s,[p.ns p.nyP p.nyF p.nb]);

        % interpolate from x(s) to get s(x)
        sav = zeros(p.ns,p.nyP,p.nyF,p.nb);
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
    
    % Distribution over (x,yP,yF,beta)
    [model.SSdist,model.sav_longgrid]...
            = find_stationary(p,model,income,prefs,xgrid.longgrid);
    
    % Cumulative distribution, sorted by saving
    temp = sortrows([model.sav_longgrid(:) model.SSdist(:)]);
    model.sav_longgrid_sort = temp(:,1);
    model.SSdist_sort = temp(:,2);
    model.SScumdist = cumsum(model.SSdist_sort);
    % unique values on cumdist and their indices (needed for interpolants)
    [model.SScumdist_unique,model.SScumdist_uniqueind] = unique(model.SScumdist,'last');
    
    % Distribution over (assets,yP_lag,yF_lag,beta_lag)
    model.asset_values = p.R * model.sav_longgrid;
    if p.Bequests == 1
        model.asset_dist = model.SSdist;
    else
        % Shift fraction of distribution to zero
        model.asset_dist = (1-p.dieprob) * model.SSdist;
        model.asset_dist(1,:,:,:) = p.dieprob * sum(model.SSdist,1);
    end
    asset_temp = sortrows([model.asset_values(:) model.asset_dist(:)]);
    model.asset_sortvalues = asset_temp(:,1);
    model.asset_dist_sort  = asset_temp(:,2);
    model.asset_cumdist    = cumsum(model.asset_dist_sort);
    % Get unique values from cumdist for interpolant
    [model.asset_cumdist_unique,model.asset_uniqueind] = unique(model.asset_cumdist,'last');
    
    % mean saving, mean assets
    model.mean_s = model.sav_longgrid(:)' * model.SSdist(:);
    model.mean_a = model.asset_values(:)' * model.asset_dist(:);
 
    % Collapse the asset distribution from (a,yP_lag,yF_lag,beta_lag) to (a,beta_lag) for norisk
    % model, and from (x,yP,yF,beta) to (x,beta)
    if p.nyP>1 && p.nyF>1
        % a
        model.assetdist_noincrisk =  sum(sum(model.asset_dist,3),2);
        % x
        model.SSdist_noincrisk    = sum(sum(model.SSdist,3),2);
    elseif (p.nyP>1 && p.nyF==1) || (p.nyP==1 && p.nyF>1)
        model.assetdist_noincrisk =  sum(model.asset_dist,2);
        model.SSdist_noincrisk    = sum(model.SSdist,2);
    elseif p.nyP==1 && p.nyF==1
        model.assetdist_noincrisk = model.asset_dist;
        model.SSdist_noincrisk    = model.SSdist;
    end
    
    % Policy functions on longgrid
    model.con_longgrid = repmat(xgrid.longgrid(:),p.nb,1)...
        - model.sav_longgrid(:) - p.savtax*max(model.sav_longgrid(:)-p.savtaxthresh,0);
    model.con_longgrid = reshape(model.con_longgrid,[p.nxlong,p.nyP,p.nyF,p.nb]);
    
    % Mean consumption
    model.mean_c = model.con_longgrid(:)' * model.SSdist(:);

    if p.Display == 1
        fprintf(' A/Y = %2.3f\n',model.mean_a/(income.meany1*p.freq));
    end
    AYdiff = model.mean_a/(income.meany1*p.freq) -  p.targetAY;

end