function [AYdiff,model] = solve_EGP_EZ(beta,p,grids,gridsKFE,prefs,income)
    
    %% CONSTRUCT EXPECTATIONS MATRIX                                     
    betagrid = beta + prefs.betagrid0;
    
    if p.IterateBeta == 1 && p.Display == 1
        msg = sprintf(' %3.3f',betagrid);
        disp([' Trying betagrid =' msg])
    end
    
    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    con = (1/p.R) * repmat(grids.x.matrix(:),p.nb,1);
    
    % initial guess for value function
    V = con;
    
    % discount factor matrix, 
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if numel(p.invies) > 1
        betastacked = speye(p.nyP*p.nyF*p.nx*p.nb) * betagrid;
    else
        betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));
        betastacked = sparse(diag(betastacked));
    end

    % Expectations operator (conditional on yT)
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if numel(p.invies) == 1
        Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));
    else
        Emat = kron(prefs.IEStrans,kron(income.ytrans,speye(p.nx)));
        invies_col = kron(p.invies',ones(p.nx*p.nyP*p.nyF,1));
        invies_col_yT = repmat(invies_col,1,p.nyT);
    end
    
    %% EGP Iteration
    iter = 1;
    cdiff = 1;
    while iter<p.max_iter && cdiff>p.tol_iter
        if iter==1
            conlast = con;
            Vlast   = V;
        else
            conlast = conupdate;
            Vlast   = Vupdate;
        end
        iter = iter + 1;
        
        % interpolate to get c(x') using c(x)
        
        % c(x) and V(x)
        conlast = reshape(conlast,[p.nx p.nyP p.nyF p.nb]);
        Vlast   = reshape(Vlast,[p.nx p.nyP p.nyF p.nb]);
        % c(x') and V(x')
        c_xp = zeros(p.nx,p.nyP,p.nyF,p.nb,p.nyT);
        V_xp = zeros(p.nx,p.nyP,p.nyF,p.nb,p.nyT);
        
        % x'(s)
        temp_sav = repmat(grids.s.matrix(:),p.nb,p.nyT);
        temp_inc = repmat(kron(income.netymat,ones(p.nx,1)),p.nb,1);
        xp_s = (1+p.r)*temp_sav + temp_inc;
        xp_s = reshape(xp_s,[p.nx p.nyP p.nyF p.nb p.nyT]);
        
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);
            coninterp = griddedInterpolant(grids.x.matrix(:,iyP,iyF),conlast(:,iyP,iyF,ib),'linear');
            c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
            Vinterp{iyP,iyF,ib} = griddedInterpolant(grids.x.matrix(:,iyP,iyF),Vlast(:,iyP,iyF,ib),'linear');
            V_xp(:,iyP,iyF,ib,:) = reshape(Vinterp{iyP,iyF,ib}(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
        end
        end
        end
        
         % reshape to take expecation over yT first
        c_xp = reshape(c_xp,[],p.nyT);
        V_xp = reshape(V_xp,[],p.nyT);

        % matrix of next period muc, muc(x',yP',yF)
        if numel(p.invies) == 1
            mucnext = c_xp.^(-p.invies) .* V_xp.^(p.invies-p.risk_aver);
        else
            mucnext = c_xp.^(-invies_col_yT) .* V_xp.^(invies_col_yT-p.risk_aver);
        end
        
        % expected muc
        savtaxrate  = (1+p.savtax.*(repmat(grids.s.matrix(:),p.nb,1)>=p.savtaxthresh));
        emuc = (1+p.r)*betastacked*Emat*mucnext*income.yTdist ./ savtaxrate;
        if p.risk_aver == 1
            ezvalnext = exp(Emat * log(V_xp) * income.yTdist);
        else
            ezvalnext = (Emat * V_xp.^(1-p.risk_aver) * income.yTdist).^(1/(1-p.risk_aver));
        end
        
        if numel(p.invies) == 1
            muc_s = emuc .* ezvalnext .^(p.risk_aver-p.invies);
            con_s = muc_s .^ (-1/p.invies);
        else
            muc_s = emuc .* ezvalnext .^(p.risk_aver-invies_col);
            con_s = muc_s .^ (-1./invies_col);
        end
        
        x_s = con_s + repmat(grids.s.matrix(:),p.nb,1)...
                        + p.savtax * max(repmat(grids.s.matrix(:),p.nb,1)-p.savtaxthresh,0);
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
        sav = max(sav,p.borrow_lim);
        xp = p.R * repmat(sav(:),1,p.nyT) ... 
                    + repmat(kron(income.netymat,ones(p.nx,1)),p.nb,1);
        xp = reshape(xp,[p.nx p.nyP p.nyF p.nb p.nyT]);
        
        conupdate = repmat(grids.x.matrix,[1 1 1 p.nb]) - sav - p.savtax * max(sav-p.savtaxthresh,0);
        
        % interpolate adjusted expected value function on x grid
        ezval_integrand = zeros(p.nx,p.nyP,p.nyF,p.nb,p.nyT);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            xp_iyP_iyF_ib = xp(:,iyP,iyF,ib,:);
            temp_iyP_iyF_ib = Vinterp{iyP,iyF,ib}(xp_iyP_iyF_ib(:)) .^ (1-p.risk_aver);
            ezval_integrand(:,iyP,iyF,ib,:) = reshape(temp_iyP_iyF_ib,[p.nx 1 1 1 p.nyT]);
        end
        end
        end
        % Take expectation over yT
        ezval_integrand = reshape(ezval_integrand,[],p.nyT) * income.yTdist;
        % Take expectation over (yP,yF,beta)
        ezval = Emat * ezval_integrand;
        if p.risk_aver == 1
            ezval = exp(ezval);
        else
            ezval = ezval .^ (1/(1-p.risk_aver));
        end

        % update value function
        ezval = reshape(ezval,p.nx,p.nyP,p.nyF,p.nb);
        Vupdate = zeros(p.nx,p.nyP,p.nyF,p.nb);
        for ib = 1:p.nb
            if numel(p.invies) == 1
                if p.invies == 1
                    Vupdate(:,:,:,ib) = conupdate(:,:,:,ib) .^ (1-betagrid(ib)) .* ezval(:,:,:,ib) .^ betagrid(ib);
                else
                	Vupdate(:,:,:,ib) = (1-betagrid(ib)) * conupdate(:,:,:,ib) .^ (1-p.invies) ...
                                    + betagrid(ib) * ezval(:,:,:,ib) .^ (1-p.invies);
                    Vupdate(:,:,:,ib) = Vupdate(:,:,:,ib) .^ (1/(1-p.invies));
                end
            else
                if p.invies(ib) == 1
                    Vupdate(:,:,:,ib) = conupdate(:,:,:,ib) .^ (1-betagrid) .* ezval(:,:,:,ib) .^ betagrid;
                else
                    invies_mat = repmat(p.invies(ib),[p.nx p.nyP p.nyF 1]);
                    Vupdate(:,:,:,ib) = (1-betagrid) * conupdate(:,:,:,ib) .^ (1-invies_mat) ...
                                    + betagrid * ezval(:,:,:,ib) .^ (1-invies_mat);
                    Vupdate(:,:,:,ib) = Vupdate(:,:,:,ib) .^ (1./(1-invies_mat));
                end
            end
        end
        
        cdiff = max(abs(conupdate(:)-conlast(:)));
        if p.Display >=1 && mod(iter,50) ==0
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
            griddedInterpolant(grids.x.matrix(:,iyP,iyF),model.sav(:,iyP,iyF,ib),'linear');
        model.coninterp{iyP,iyF,ib} = ...
            griddedInterpolant(grids.x.matrix(:,iyP,iyF),model.con(:,iyP,iyF,ib),'linear');    
    end
    end
    end
    
    %% DISTRIBUTION
     
    model = find_stationary_adist(p,model,income,prefs,gridsKFE);
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
        model.xdist_noincrisk = sum(sum(model.xdist,3),2);
    elseif (p.nyP>1 && p.nyF==1) || (p.nyP==1 && p.nyF>1)
        model.adist_noincrisk =  sum(model.adist,2);
        model.xdist_noincrisk = sum(model.xdist,2);
    elseif p.nyP==1 && p.nyF==1
        model.adist_noincrisk = model.adist;
        model.xdist_noincrisk = model.xdist;
    end

    % consumption policy function on xgrid
    model.con_x = model.xvals - model.sav_x - p.savtax*max(model.sav_x-p.savtaxthresh,0);
           
    % mean saving, mean assets
    model.mean_a = model.adist(:)' * grids.a.matrix(:);
 
    
    if p.Display == 1
        fprintf(' A/Y = %2.3f\n',model.mean_a/(income.meany1*p.freq));
    end
    AYdiff = model.mean_a/(income.meany1*p.freq) -  p.targetAY;
    
end