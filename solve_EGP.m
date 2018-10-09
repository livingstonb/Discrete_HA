function [AYdiff,model] = solve_EGP(beta,p,xgrid,sgrid,prefs,income)
    % This function performs the method of endogenous grid points to find
    % saving and consumption policy functions. It also calls 
    % find_stationary() to compute the stationary distribution over states 
    % via direct methods (rather than simulations) and stores the results 
    % in the 'model' structure.
    
    %% CONSTRUCT EXPECTATIONS MATRIX                                     
    betagrid = beta + prefs.betagrid0;
    
    if p.IterateBeta == 1
        msg = sprintf(' %3.3f',betagrid);
        disp([' Trying betagrid =' msg])
    end

    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    con = p.r * repmat(xgrid.orig_wide(:),p.nb,1);

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
        conlast_wide = reshape(conlast,[p.ns p.nyP p.nyF p.nb]);
        % c(x')
        c_xp = zeros(p.ns,p.nyP,p.nyF,p.nb,p.nyT);
        
        % x'(s)
        temp_sav_wide = repmat(sgrid.wide(:),p.nb,p.nyT);
        temp_inc_wide = repmat(kron(income.netymat,ones(p.ns,1)),p.nb,1);
        xp_s_wide = (1+p.r)*temp_sav_wide + temp_inc_wide;
        xp_s_wide = reshape(xp_s_wide,[p.ns p.nyP p.nyF p.nb p.nyT]);

        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            coninterp = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF),conlast_wide(:,iyP,iyF,ib),'linear');
            xp_s_wide_ib_iyF_iyP = xp_s_wide(:,iyP,iyF,ib,:);
            c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(xp_s_wide_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
        end
        end
        end
        
        % reshape to take expecation over yT first
        c_xp        = reshape(c_xp,[],p.nyT);
        xp_s_wide   = reshape(xp_s_wide,[],p.nyT);

        % matrix of next period muc, muc(x',yP',yF)
        mucnext     = prefs.u1(c_xp) - p.temptation/(1+p.temptation) * prefs.u1(xp_s_wide);
        
        % now find muc this period as a function of s:
        % variables defined for each (x,yP,yF,beta) in state space,
        % column vecs of length p.nx*p.nyP*p.nyF*p.nb
        muc_savtaxrate  = (1+p.savtax.*(repmat(sgrid.wide(:),p.nb,1)>=p.savtaxthresh));
        muc_consumption = (1+p.r)*betastacked*Emat*(mucnext*income.yTdist);
        muc_bequest     = p.dieprob*prefs.beq1(repmat(sgrid.wide(:),p.nb,1));
        
        % muc(s(x,yP,yF,beta))
        muc_s           = (1-p.dieprob) * muc_consumption ./ muc_savtaxrate...
                                                + p.dieprob * muc_bequest;
                
        % c(s)
        con_s       = prefs.u1inv(muc_s);
        % x(s) = s + stax + c(s)
        x_s_wide    = repmat(sgrid.wide(:),p.nb,1)...
                        + p.savtax * max(repmat(sgrid.wide(:),p.nb,1)-p.savtaxthresh,0)...
                        + con_s;
        x_s_wide = reshape(x_s_wide,[p.ns p.nyP p.nyF p.nb]);

        % interpolate from x(s) to get s(x)
        sav_wide = zeros(p.ns,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            savinterp = griddedInterpolant(x_s_wide(:,iyP,iyF,ib),sgrid.wide(:,iyP,iyF),'linear');
            sav_wide(:,iyP,iyF,ib) = savinterp(xgrid.orig_wide(:,iyP,iyF)); 
        end
        end
        end

        % deal with borrowing limit
        sav_wide(sav_wide<p.borrow_lim) = p.borrow_lim;

        % updated consumption function, column vec length of
        % length p.nx*p.nyP*p.nyF*p.nb
        conupdate = repmat(xgrid.orig_wide(:),p.nb,1) - sav_wide(:)...
                            - p.savtax * max(sav_wide(:)-p.savtaxthresh,0);

        cdiff = max(abs(conupdate-conlast));
        if mod(iter,50) ==0
            disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end

    end
    
    if cdiff>p.tol_iter
        % EGP did not converge, don't find stationary distribution
        AYdiff = 100000;
        model.EGP_cdiff = cdiff;
        return
    end

    model.sav_wide = sav_wide;
    model.sav = sav_wide(:);
    model.con = conupdate;
    model.con_wide = reshape(conupdate,[p.nx p.nyP p.nyF p.nb]);
    model.EGP_cdiff = cdiff;
    
    % create interpolants from optimal policy functions
    model.savinterp = cell(p.nyP,p.nyF,p.nb);
    model.coninterp = cell(p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        model.savinterp{iyP,iyF,ib} = ...
            griddedInterpolant(xgrid.orig_wide(:,iyP,iyF),model.sav_wide(:,iyP,iyF,ib),'linear');
        model.coninterp{iyP,iyF,ib} = ...
            griddedInterpolant(xgrid.orig_wide(:,iyP,iyF),model.con_wide(:,iyP,iyF,ib),'linear');
    end
    end
    end


    %% DISTRIBUTION

    SkipStationary = 0;
    [model.SSdist,model.statetrans,model.sav_longgrid_wide]...
            = find_stationary(p,model,income,prefs,xgrid.longgrid_wide,SkipStationary);

    % SS probability of residing in each state
    model.SSdist_wide = reshape(model.SSdist,[p.nxlong,p.nyP,p.nyF,p.nb]);
    
    % SS distribution of assets
    model.asset_values = p.R * model.sav_longgrid_wide;
    if p.WealthInherited == 1
        model.asset_dist = model.SSdist_wide;
    else
        % Shift fraction of distribution to zero
        model.asset_dist = (1-p.dieprob) * model.SSdist_wide;
        model.asset_dist(1,:,:,:) = p.dieprob * sum(model.SSdist_wide,1);
    end
    asset_temp = sortrows([model.asset_values(:) model.asset_dist(:)]);
    model.asset_sortvalues = asset_temp(:,1);
    model.asset_dist_sort  = asset_temp(:,2);
    model.asset_cumdist    = cumsum(model.asset_dist_sort);
    % Get unique values from cumdist for interpolant
    [model.asset_cumdist_unique,model.asset_uniqueind] = unique(model.asset_cumdist,'last');
    
    % mean saving
    model.mean_s = model.sav_longgrid_wide(:)' * model.SSdist;
    
    % policy functions on longgrid
    model.sav_longgrid = model.sav_longgrid_wide(:);
    model.con_longgrid = repmat(xgrid.longgrid_wide(:),p.nb,1)...
        - model.sav_longgrid(:) - p.savtax*max(model.sav_longgrid(:)-p.savtaxthresh,0);
    model.con_longgrid_wide = reshape(model.con_longgrid,[p.nxlong,p.nyP,p.nyF,p.nb]);
    model.mean_c = model.con_longgrid_wide(:)' * model.SSdist;
    
    % mean assets
    if p.WealthInherited == 1
        model.mean_a = p.R * model.mean_s;
    else
        model.mean_a = (1-p.dieprob) * p.R * model.mean_s;
    end
    
    % cumulative distribution, sorted by saving
    temp = sortrows([model.sav_longgrid model.SSdist]);
    model.sav_longgrid_sort = temp(:,1);
    model.SSdist_sort = temp(:,2);
    model.SScumdist = cumsum(model.SSdist_sort);
    
    %
    
    % Collapse the distribution from (x,yP,yF,beta) to (x,beta) for norisk
    % model
    if p.nyP>1 && p.nyF>1
        model.SSdist_noincrisk =  sum(model.SSdist_wide,[2 3]);
    elseif (p.nyP>1 && p.nyF==1) || (p.nyP==1 & p.nyF>1)
        model.SSdist_noincrisk =  sum(model.SSdist_wide,2);
    elseif p.nyP==1 && p.nyF==1
        model.SSdist_noincrisk = model.SSdist_wide;
    end
    
    % unique values on cumdist and their indices (needed for interpolants)
    [model.SScumdist_unique,model.SScumdist_uniqueind] = unique(model.SScumdist,'last');

    fprintf(' A/Y = %2.3f\n',model.mean_a/(income.meany*p.freq));
    AYdiff = model.mean_a/(income.meany*p.freq) -  p.targetAY;

end