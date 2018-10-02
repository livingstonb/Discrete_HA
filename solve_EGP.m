function [AYdiff,model,xgridm] = solve_EGP(beta,p,xgrid,sgrid,prefs,...
                                            ergodic_tol,income,gridsize)
    % This function performs the method of endogenous grid points to find
    % saving and consumption policy functions. It also computes the
    % stationary distribution over states via direct methods (rather than
    % simulations) and stores the results in the 'model' structure.
    
    % 'ergodic_tol' specifies the tolerance used to find the ergodic
    % distribution from the transition matrix
    

    % 'gridsize' specifies the number of points in the asset space used to 
    % find the ergodic distribution, which may be larger than the number of
    % points used to find the policy functions
    
    %% CONSTRUCT EXPECTATIONS MATRIX                                     
    if  p.nb == 1
        betagrid = beta;
    elseif p.nb ==2 
        betagrid = [beta-p.betawidth;beta+p.betawidth];
    end

    if p.IterateBeta == 1
        disp(['Trying betagrid = ' num2str(betagrid)])
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
    model.Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));

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
        xp_s_wide = (1+p.r)*repmat(sgrid.wide(:),p.nb,p.nyT) + repmat(kron(income.netymat,ones(p.ns,1)),p.nb,1);
        xp_s_wide = reshape(xp_s_wide,[p.ns p.nyP p.nyF p.nb p.nyT]);

        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            coninterp = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF),conlast_wide(:,iyP,iyF,ib),'linear','linear');
            c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(reshape(xp_s_wide(:,iyP,iyF,ib,:),[],1)),[],1,1,1,p.nyT);
        end
        end
        end
        
        % reshape to take expecation yT first
        c_xp        = reshape(c_xp,[],p.nyT);
        xp_s_wide   = reshape(xp_s_wide,[],p.nyT);

        % matrix of next period muc, muc(x',yP',yF)
        mucnext     = prefs.u1(c_xp) - p.temptation/(1+p.temptation) * prefs.u1(xp_s_wide);
        
        % now find muc this period as a function of s:
        % variables defined for each (x,yP,yF,beta) in state space,
        % column vecs of length p.nx*p.nyP*p.nyF*p.nb
        muc_savtaxrate  = (1+p.savtax.*(repmat(sgrid.wide(:),p.nb,1)>=p.savtaxthresh));
        muc_consumption = (1+p.r)*betastacked*(model.Emat*(mucnext*income.yTdist));
        muc_bequest     = p.dieprob*prefs.beq1(repmat(sgrid.wide(:),p.nb,1));
        % muc(s(x,yP,yF,beta))
        muc_s           = (1-p.dieprob) * muc_consumption .* muc_savtaxrate...
                            + p.dieprob * muc_bequest;
                
        % c(s)
        con_s       = prefs.u1inv(muc_s);
        % x(s) = s + stax + c(s)
        x_s_wide    = repmat(sgrid.wide(:),p.nb,1)...
                        + p.savtax * max(repmat(sgrid.wide(:),p.nb,1)-p.savtaxthresh,0)...
                        + con_s;

        % interpolate from x(s) to get s(x), interpolate for each (beta,yP,yF)
        % separately
        x_s_wide = reshape(x_s_wide,[p.ns p.nyP p.nyF p.nb]);
        sav_wide = zeros(p.ns,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            savinterp = griddedInterpolant(x_s_wide(:,iyP,iyF,ib),sgrid.wide(:,iyP,iyF),'linear','linear');
            sav_wide(:,iyP,iyF,ib) = savinterp(xgrid.orig_wide(:,iyP,iyF)); 
        end
        end
        end

        % deal with borrowing limit
        sav_wide(sav_wide<p.borrow_lim) = p.borrow_lim;

        % updated consumption function, column vec length of
        % p.nx*p.nyP*p.nyF*p.nb
        conupdate = repmat(xgrid.orig_wide(:),p.nb,1) - sav_wide(:) - p.savtax * max(sav_wide(:)-p.savtaxthresh,0);

        cdiff = max(abs(conupdate-conlast));
        if mod(iter,50) ==0
            disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end

    end
    
    if cdiff>p.tol_iter
        % EGP did not converge, don't find stationary distribution
        AYdiff = 100000;
        model.EGP_cdiff = cdiff;
        xgridm = [];
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
    
    [model.SSdist,model.statetrans,model.sav_longgrid_wide]...
            = find_stationary(p,model,income,prefs,xgridm,ergodic_tol);

    % SS probability of residing in each state
    model.SSdist_wide = reshape(model.SSdist,[gridsize,p.nyP,p.nyF,p.nb]);

    % SS wealth/gross income ratio
    model.mean_s = model.sav_longgrid_wide(:)' * model.SSdist;
    model.mean_a = p.R * model.mean_s;

    % policy functions
    model.sav_longgrid      = model.sav_longgrid_wide(:);
    if p.WealthInherited == 0
        model.a_longgrid        = (1 - p.dieprob) * p.R * model.sav_longgrid;
    else
        model.a_longgrid        = p.R * model.sav_longgrid;
    end
    model.con_longgrid      = repmat(xgridm(:),p.nb,1) - model.sav_longgrid(:) - p.savtax*max(model.sav_longgrid-p.savtaxthresh,0);
    model.con_longgrid_wide = reshape(model.con_longgrid,[gridsize,p.nyP,p.nyF,p.nb]);
    
    % cumulative distribution
    temp = sortrows([model.sav_longgrid model.SSdist]);
    model.sav_longgrid_sort = temp(:,1);
    model.a_longgrid_sort = (1 - p.dieprob) * p.R * model.sav_longgrid_sort;
    model.SSdist_sort = temp(:,2);
    model.SScumdist = cumsum(model.SSdist_sort);
    
    % unique values on cumdist and their indices (needed for interpolants)
    [model.SScumdist_unique,model.SScumdist_uniqueind] = unique(model.SScumdist,'last');
    
    mean_s = model.SSdist' * model.sav_longgrid;

    fprintf(' A/Y = %2.3f\n',model.mean_a/income.meany);
    %AYdiffsq = (mean_s/meany - targetAY)^2;
    AYdiff = model.mean_a/income.meany -  p.targetAY;

end