function [AYdiff,model] = solve_EGP_EZ(beta,p,grids,gridsKFE,prefs,income)
    
    %% CONSTRUCT EXPECTATIONS MATRIX                                     
    betagrid = beta + prefs.betagrid0;
    
    if p.IterateBeta == 1 && p.Display == 1
        msg = sprintf(' %3.3f',betagrid);
        disp([' Trying betagrid =' msg])
    end
    
    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    if p.r < 0.001
        extra = 0.002;
    else
        extra = 0;
    end
    con = (p.r + extra) * repmat(grids.x.matrix(:),p.nb,1);
    
    % initial guess for value function
    V = con;
    
    % discount factor matrix, 
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if (numel(p.invies) > 1) || (numel(p.risk_aver) > 1)
        betastacked = speye(p.nyP*p.nyF*p.nx*p.nb) * betagrid;
    else
        betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));
        betastacked = sparse(diag(betastacked));
    end

    % Expectations operator (conditional on yT)
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb   
    if numel(p.invies) > 1
        Emat = kron(prefs.ztrans,kron(income.ytrans,speye(p.nx)));
        invies_col = kron(p.invies',ones(p.nx*p.nyP*p.nyF,1));
        risk_aver_col = p.risk_aver;
        invies_col_yT = repmat(invies_col,1,p.nyT);
    elseif numel(p.risk_aver) > 1
        Emat = kron(prefs.ztrans,kron(income.ytrans,speye(p.nx)));
        risk_aver_col = kron(p.risk_aver',ones(p.nx*p.nyP*p.nyF,1));
        invies_col = p.invies;
        risk_aver_col_yT = repmat(risk_aver_col,1,p.nyT);
    else
        Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));
        risk_aver_col = p.risk_aver;
        invies_col = p.invies;
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
        if numel(p.invies) > 1
            mucnext = c_xp.^(-invies_col_yT) .* V_xp.^(invies_col_yT-p.risk_aver);
        elseif numel(p.risk_aver) > 1
            mucnext = c_xp.^(-p.invies) .* V_xp.^(p.invies-risk_aver_col_yT);
        else
            mucnext = c_xp.^(-p.invies) .* V_xp.^(p.invies-p.risk_aver);
        end
        
        % expected muc
        savtaxrate  = (1+p.savtax.*(repmat(grids.s.matrix(:),p.nb,1)>=p.savtaxthresh));
        mu_cons = (1+p.r)*betastacked*Emat*mucnext*income.yTdist ./ savtaxrate;
        mu_bequest = prefs.beq1(repmat(grids.s.matrix(:),p.nb,1));
        emuc = (1-p.dieprob) * mu_cons + p.dieprob * mu_bequest;
        
        if numel(p.risk_aver) > 1
            ezvalnext_ra_equal1 = exp(Emat * log(V_xp) * income.yTdist)...
                .^ (risk_aver_col-invies_col);

            ezvalnext_ra_nequal1 = (Emat * V_xp.^(1-risk_aver_col) * income.yTdist)...
                .^ ((risk_aver_col-invies_col)./(1-risk_aver_col));

            ezvalnext = zeros(p.nx*p.nyP*p.nyF*p.nb,1);
            ezvalnext(risk_aver_col==1) = ezvalnext_ra_equal1(risk_aver_col==1);
            ezvalnext(risk_aver_col~=1) = ezvalnext_ra_nequal1(risk_aver_col~=1);
        else
            if p.risk_aver == 1
                ezvalnext = exp(Emat * log(V_xp) * income.yTdist)...
                    .^ (risk_aver_col - invies_col);
            else
                ezvalnext = (Emat * V_xp .^ (1-risk_aver_col) * income.yTdist)...
                    .^ ((risk_aver_col - invies_col)./(1 - risk_aver_col));
            end
        end
        
        muc_s = emuc .* ezvalnext;
        con_s = muc_s .^ (-1./invies_col);
        
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

        index_to_extend = 1*(p.nyF==1) + 2*(p.nyF>1);
        xp = p.R * repmat(sav(:),1,p.nyT) ... 
                    + repmat(kron(income.netymat,ones(p.nx,1)),p.nb,index_to_extend);
        xp = reshape(xp,[p.nx p.nyP p.nyF p.nb p.nyT]);
        
        conupdate = repmat(grids.x.matrix,[1 1 1 p.nb]) - sav - p.savtax * max(sav-p.savtaxthresh,0);
        
        % interpolate adjusted expected value function on x grid
        ezval_integrand = zeros(p.nx,p.nyP,p.nyF,p.nb,p.nyT);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            xp_iyP_iyF_ib = xp(:,iyP,iyF,ib,:);
            temp_iyP_iyF_ib = Vinterp{iyP,iyF,ib}(xp_iyP_iyF_ib(:)) .^ (1-p.risk_aver(ib));
            ezval_integrand(:,iyP,iyF,ib,:) = reshape(temp_iyP_iyF_ib,[p.nx 1 1 1 p.nyT]);
        end
        end
        end
        % Take expectation over yT
        ezval_integrand = reshape(ezval_integrand,[],p.nyT) * income.yTdist;
        % Take expectation over (yP,yF,beta)
        ezval = Emat * ezval_integrand;

        ezval_ra_equal1 = (risk_aver_col==1) .* exp(ezval);
        ezval_ra_nequal1 = (risk_aver_col~=1) .* ezval .^ (1./(1-risk_aver_col));
        
        ezval = zeros(p.nx*p.nyP*p.nyF*p.nb,1);
        ezval(risk_aver_col==1) = ezval_ra_equal1(risk_aver_col==1);
        ezval(risk_aver_col~=1) = ezval_ra_nequal1(risk_aver_col~=1);

        % update value function
        ezval = reshape(ezval,p.nx,p.nyP,p.nyF,p.nb);
        Vupdate = zeros(p.nx,p.nyP,p.nyF,p.nb);
        for ib = 1:p.nb
            if numel(p.invies) == 1

                if numel(p.risk_aver) == 1
                    ibeta = ib; % possible beta heterogeneity
                else
                    ibeta = 1;
                end

                if p.invies == 1
                    Vupdate(:,:,:,ib) = conupdate(:,:,:,ib) .^ (1-betagrid(ibeta)) .* ezval(:,:,:,ib) .^ betagrid(ibeta);
                else
                	Vupdate(:,:,:,ib) = (1-betagrid(ibeta)) * conupdate(:,:,:,ib) .^ (1-p.invies) ...
                                    + betagrid(ibeta) * ezval(:,:,:,ib) .^ (1-p.invies);
                    Vupdate(:,:,:,ib) = Vupdate(:,:,:,ib) .^ (1/(1-p.invies));
                end
            elseif numel(p.invies) > 1
                if p.invies(ib) == 1
                    Vupdate(:,:,:,ib) = conupdate(:,:,:,ib) .^ (1-betagrid) .* ezval(:,:,:,ib) .^ betagrid;
                else
                    Vupdate(:,:,:,ib) = (1-betagrid) * conupdate(:,:,:,ib) .^ (1-p.invies(ib)) ...
                                    + betagrid * ezval(:,:,:,ib) .^ (1-p.invies(ib));
                    Vupdate(:,:,:,ib) = Vupdate(:,:,:,ib) .^ (1/(1-p.invies(ib)));
                end
            end
        end
        
        assert(sum(~isfinite(Vupdate(:)))==0)
        
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
    model.V = reshape(Vupdate,[p.nx p.nyP p.nyF p.nb]);
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
    model.mean_a = model.adist(:)' * gridsKFE.a.matrix(:);
 
    
    if p.Display == 1
        fprintf(' A/Y = %2.3f\n',model.mean_a/(income.meany1*p.freq));
    end
    AYdiff = model.mean_a/(income.meany1*p.freq) -  p.targetAY;
    
end