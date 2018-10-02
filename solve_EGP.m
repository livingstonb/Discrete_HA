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
    fprintf(' Computing state-to-state transition probabilities... \n');
    savm = model.sav_wide;

    % Create (potentially) longer grid
    xgridm = linspace(0,1,gridsize)';
    xgridm = repmat(xgridm,[1 p.nyP p.nyF]) .^(1/p.xgrid_par);
    % for each (yP,yF), add nety income for smallest yT value to grid
    min_netymat = reshape(income.netymat(:,1),[1 p.nyP p.nyF]);
    min_netymat = repmat(min_netymat,[gridsize 1 1]);
    xgridm = p.borrow_lim + min_netymat + (p.xmax-p.borrow_lim)*xgridm;

    % Interpolate policy functions onto new grid
    savlong = zeros(gridsize,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        savlong(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(xgridm(:,iyP,iyF));
    end
    end
    end
    savm = savlong;
    NN = gridsize * p.nyP * p.nyF * p.nb;
    nn = gridsize;

    netymatm = reshape(income.netymat,[p.nyP p.nyF p.nyT]);
    
    trans = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    grid_probabilities = zeros(NN,NN);
    col = 1;
    for ib2 = 1:p.nb
    for iyF2 = 1:p.nyF
    for iyP2 = 1:p.nyP
        fspace = fundef({'spli',xgridm(:,iyP2,iyF2),0,1});
        % xprime if no death
        xp_live = (1+p.r)*repmat(savm(:),p.nyT,1) + ...
            kron(squeeze(netymatm(iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));

        % xprime if death (i.e. saving in past period set to 0)
        if p.WealthInherited == 0
            xp_death = kron(squeeze(netymatm(iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
        else
            xp_death = xp_live;
        end
        
        % set probabilities equal to 1 at grid endpt when xp is off the grid
        idx_xpl_max = xp_live>=max(xgridm(:,iyP2,iyF2));
        idx_xpl_min = xp_live<=min(xgridm(:,iyP2,iyF2));
        
        idx_xpd_max = xp_death>=max(xgridm(:,iyP2,iyF2));
        idx_xpd_min = xp_death<=min(xgridm(:,iyP2,iyF2));
        
        interpl = funbas(fspace,xp_live);
        interpl(idx_xpl_max | idx_xpl_min,:) = 0;
        interpl(idx_xpl_max,end) = 1;
        interpl(idx_xpl_min,1) = 1;
        
        interpd = funbas(fspace,xp_death);
        interpd(idx_xpd_max | idx_xpd_min,:) = 0;
        interpd(idx_xpd_max,end) = 1;
        interpd(idx_xpd_min,1) = 1;
        
        interp = (1-p.dieprob) * interpl + p.dieprob * interpd;
        
        % if not setting probabilites to 1 at grid endpts
        % interp = (1-p.dieprob) * funbas(fspace,xp_live) + p.dieprob * funbas(fspace,xp_death);
        interp = reshape(interp,[],p.nyT*nn);
        % Multiply by yT distribution
        newcolumn = interp * kron(speye(nn),income.yTdist);
        % Multiply by (beta,yF,yP) distribution
        newcolumn = bsxfun(@times,kron(trans(:,col),ones(nn,1)),newcolumn);

        grid_probabilities(:,nn*(col-1)+1:nn*col) = newcolumn;
        col = col + 1;
    end
    end
    end
    
    model.statetrans = grid_probabilities;

    % SS probability of residing in each state
    fprintf(' Finding ergodic distribution...\n');
    model.SSdist = full(ergodicdist(sparse(grid_probabilities),1,ergodic_tol));
    model.SSdist_wide = reshape(model.SSdist,[nn,p.nyP,p.nyF,p.nb]);

    % SS wealth/gross income ratio
    model.mean_s = savm(:)' * model.SSdist;
    model.mean_a = p.R * model.mean_s;

    % policy functions
    model.sav_longgrid      = savm(:);
    model.sav_longgrid_wide = reshape(model.sav_longgrid,[nn,p.nyP,p.nyF,p.nb]);
    if p.WealthInherited == 0
        model.a_longgrid        = (1 - p.dieprob) * p.R * model.sav_longgrid;
    else
        model.a_longgrid        = p.R * model.sav_longgrid;
    end
    model.con_longgrid      = repmat(xgridm(:),p.nb,1) - savm(:) - p.savtax*max(savm(:)-p.savtaxthresh,0);
    model.con_longgrid_wide = reshape(model.con_longgrid,[nn,p.nyP,p.nyF,p.nb]);
    
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