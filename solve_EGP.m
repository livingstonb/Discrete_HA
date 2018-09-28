function [AYdiff,model,xgridm] = solve_EGP(beta,p,xgrid,sgrid,prefs,...
                                            ergodic_tol,income,gridsize)


    if  p.nb == 1
        betagrid = beta;
    elseif p.nb ==2 
        betagrid = [beta-p.betawidth;beta+p.betawidth];
    end

    if p.IterateBeta == 1
        disp(['Trying betagrid = ' num2str(betagrid)])
    end

    % initial guess for consumption function
    con = p.r * xgrid.orig_wide(:);

    % discount factor matrix
    betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));

    % Expectations operator (conditional on yT)
    Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));

    iter = 1;
    cdiff = 1;
    %% EGP ITERATION
    while iter<p.max_iter && cdiff>p.tol_iter
        if iter==1
            conlast = con;
        else
            conlast = conupdate;
        end
        iter = iter + 1;

        % interpolate to get c(x') using c(x)
        conlast_wide = reshape(conlast,[p.ns p.nyP p.nyF p.nb]);

        x_s_wide = (1+p.r)*repmat(sgrid.wide(:),1,p.nyT) + kron(income.netymat,ones(p.ns,1));
        x_s_wide = reshape(x_s_wide,[p.ns p.nyP p.nyF p.nb p.nyT]);

        % c(x')
        c_xp = zeros(p.ns,p.nyP,p.nyF,p.nb,p.nyT);
        for ib = 1:p.nb
        for iyF  = 1:p.nyF
        for iyP = 1:p.nyP
            coninterp = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),conlast_wide(:,iyP,iyF,ib),'linear','linear');
            c_xp_temp = coninterp(reshape(x_s_wide(:,iyP,iyF,ib,:),[],1));
            c_xp(:,iyP,iyF,ib,:) = reshape(c_xp_temp,[],1,1,1,p.nyT);
        end
        end
        end
        c_xp = reshape(c_xp,[],p.nyT);
        x_s_wide = reshape(x_s_wide,[],p.nyT);

        mucnext  = prefs.u1(c_xp) - p.temptation/(1+p.temptation)*prefs.u1(x_s_wide);
        % muc this period as a function of s
        muc_s = (1-p.dieprob)*(1+p.r)*betastacked.*(Emat*(mucnext*income.yTdist))./(1+p.savtax.*(sgrid.wide(:)>=p.savtaxthresh))...
            + p.dieprob*prefs.beq1(sgrid.wide(:));
        % _wide variables have dimension nx by nyP*nyF*nb, or nx by N/nx

        % consumption as a function of s
        con_s = prefs.u1inv(muc_s);
        % cash-in-hand (x) as a function of s
        x_s_wide = sgrid.wide(:) + p.savtax * max(sgrid.wide(:)-p.savtaxthresh,0) + con_s;

        % interpolate from x(s) to get s(x), interpolate for each (beta,yP,yF)
        % separately
        x_s_wide = reshape(x_s_wide,[p.ns p.nyP p.nyF p.nb]);
        sav_wide = zeros(p.ns,p.nyP,p.nyF,p.nb);
        sgrid_reshaped = reshape(sgrid.wide,[p.ns p.nyP p.nyF p.nb]);
        for ib = 1:p.nb
        for iyF  = 1:p.nyF
        for iyP = 1:p.nyP
            savinterp = griddedInterpolant(x_s_wide(:,iyP,iyF,ib),sgrid_reshaped(:,iyP,iyF,ib),'linear','linear');
            sav_wide(:,iyP,iyF,ib) = savinterp(xgrid.orig_wide(:,iyP,iyF,ib)); 
        end
        end
        end

        % deal with borrowing limit
        sav_wide(sav_wide<p.borrow_lim) = p.borrow_lim;

        conupdate = xgrid.orig_wide(:) - sav_wide(:) - p.savtax * max(sav_wide(:)-p.savtaxthresh,0);

        cdiff = max(abs(conupdate-conlast));
        if mod(iter,50) ==0
            disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end

    end

    model.sav_wide = sav_wide;
    model.sav = sav_wide(:);
    model.con = conupdate;
    model.con_wide = reshape(conupdate,[p.nx p.nyP p.nyF p.nb]);
    model.EGP_cdiff = cdiff;


    %% DISTRIBUTION
    fprintf(' Computing state-to-state transition probabilities... \n');
    savm = model.sav_wide;

    % Create grid
    xgridm = linspace(0,1,gridsize)';
    xgridm = repmat(xgridm,[1 p.nyP p.nyF]) .^(1/p.xgrid_par);
    income.netymatm = reshape(income.netymat,[1 p.nyP p.nyF p.nyT]);
    income.netymatm = repmat(income.netymatm,[gridsize 1 1 1]);
    xgridm = p.borrow_lim + min(income.netymatm,[],4) + (p.xmax-p.borrow_lim)*xgridm;

    savlong = zeros(gridsize,p.nyP,p.nyF,p.nb);
    savinterp = cell(p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        savinterp{iyP,iyF,ib} = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),savm(:,iyP,iyF,ib),'linear','linear');
        savlong(:,iyP,iyF,ib) = savinterp{iyP,iyF,ib}(xgridm(:,iyP,iyF));
    end
    end
    end
    savm = savlong;
    NN = gridsize * p.nyP * p.nyF * p.nb;
    nn = gridsize;

    trans = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    grid_probabilities = zeros(NN,NN);
    col = 1;
    for ib2 = 1:p.nb
    for iyF2 = 1:p.nyF
    for iyP2 = 1:p.nyP
        fspace = fundef({'spli',xgridm(:,iyP2,iyF2,ib2),0,1});
        % xprime if no death
        xp_live = (1+p.r)*repmat(savm(:),p.nyT,1) + ...
            kron(squeeze(income.netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
        xp_live = max(xp_live,min(xgridm(:,iyP2,iyF2,ib2)));
        xp_live = min(xp_live,p.xmax);

        % xprime if death (i.e. saving in past period set to 0)
        if p.WealthInherited == 0
            xp_death = kron(squeeze(income.netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
        else
            xp_death = xp_live;
        end
        xp_death = max(xp_death,min(xgridm(:,iyP2,iyF2,ib2)));
        xp_death = min(xp_death,p.xmax);
        
% set probabilities equal to 1 at grid endpt where xp is off the grid
%         idx_xpl_max = xp_live==p.xmax;
%         idx_xpl_min = xp_live==p.borrow_lim;
%         
%         idx_xpd_max = xp_death==p.xmax;
%         idx_xpd_min = xp_death==p.borrow_lim;
%         
%         interp = funbas(fspace,xp_live);
%         interp(idx_xpl_max,:) = 0;
%         interp(idx_xpl_max,end) = 1;
%         interp(idx_xpl_min,:) = 0;
%         interp(idx_xpl_min,1) = 1;
%         
%         interpd = funbas(fspace,xp_death);
%         interpd(idx_xpd_max,:) = 0;
%         interpd(idx_xpd_max,end) = 1;
%         interpd(idx_xpd_min,:) = 0;
%         interpd(idx_xpd_min,1) = 1;
        
%         interp = (1-p.dieprob) * interp + p.dieprob * interpd;
        
        interp = (1-p.dieprob) * funbas(fspace,xp_live) + p.dieprob * funbas(fspace,xp_death);
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

    % SS probability of residing in each state
    fprintf(' Finding ergodic distribution...\n');
    model.SSdist = full(ergodicdist(sparse(grid_probabilities),1,ergodic_tol));
    model.SSdist_wide = reshape(model.SSdist,[nn,p.nyP,p.nyF,p.nb]);

    % SS wealth/gross income ratio
    model.mean_s = savm(:)' * model.SSdist;

    % policy functions
    model.sav_longgrid = savm(:);
    model.sav_longgrid_wide = reshape(model.sav_longgrid,[nn,p.nyP,p.nyF,p.nb]);
    model.con_longgrid = xgridm(:) - savm(:) - p.savtax*max(savm(:)-p.savtaxthresh,0);
    model.con_longgrid_wide = reshape(model.con_longgrid,[nn,p.nyP,p.nyF,p.nb]);
    
    % interpolants
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        model.savinterp{iyP,iyF,ib} = ...
            griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),model.sav_wide(:,iyP,iyF,ib),'linear');
        model.coninterp{iyP,iyF,ib} = ...
            griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),model.con_wide(:,iyP,iyF,ib),'linear');
    end
    end
    end
    
    % cumulative distribution
    temp = sortrows([model.sav_longgrid model.SSdist]);
    model.sav_longgrid_sort = temp(:,1);
    model.SSdist_sort = temp(:,2);
    model.SScumdist = cumsum(model.SSdist_sort);
    
    mean_s = model.SSdist' * model.sav_longgrid;

    fprintf(' A/Y = %2.3f\n',mean_s/income.meany);
    %AYdiffsq = (mean_s/meany - targetAY)^2;
    AYdiff = mean_s/income.meany -  p.targetAY;

end