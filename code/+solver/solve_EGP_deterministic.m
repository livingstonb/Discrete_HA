function norisk = solve_EGP_deterministic(p,grids,heterogeneity,income,direct_results)
    % This function uses the method of endogenous grid points to find the
    % policy functions of the deterministic model. Output is in the
    % 'norisk' structure.

    coninterp = cell(1,p.nb);
    sav = zeros(p.nx,p.nb);
    mucnext= zeros(p.nx,p.nb);
    
    betagrid = direct_results.beta + heterogeneity.betagrid0;
    
    if numel(p.r) > 1
        Emat = kron(heterogeneity.rtrans,kron(income.ytrans,speye(p.nx)));
        r_col = kron(p.r',ones(p.nx,1));
        r_mat = reshape(r_col,[p.nx,numel(p.r)]);
    else
        r_mat = p.r;
    end

    if numel(p.risk_aver) > 1
        risk_aver_mat = kron(p.risk_aver,ones(p.nx,1));
    end

    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    if p.temptation > 0.005
        extra = 0.5;
    else
        extra = 0;
    end
    
    con = (r_mat .* (r_mat>=0.001) + 0.001 * (r_mat<0.001) + extra) .* repmat(grids.s.vec,1,p.nb) + income.meannety1;

    iter = 0;
    cdiff = 1000;
    while iter <= p.max_iter && cdiff>p.tol_iter
        iter = iter + 1;
        
        conlast = con;
        
        for ib = 1:p.nb
            coninterp{ib} = griddedInterpolant(grids.x.vec_norisk,conlast(:,ib),'linear');
            
            % cash-on-hand is just Rs + meany
            if numel(p.r) > 1
                mucnext(:,ib) = aux.utility1(p.risk_aver,coninterp{ib}(p.R(ib)*grids.s.vec + income.meannety1))...
                                - p.temptation/(1+p.temptation) * aux.utility1(p.risk_aver,p.R(ib)*grids.s.vec + income.meannety1);
            elseif numel(p.risk_aver) > 1
                mucnext(:,ib) = aux.utility1(risk_aver_mat(:,ib),coninterp{ib}(p.R.*grids.s.vec + income.meannety1))...
                                - p.temptation/(1+p.temptation) * aux.utility1(risk_aver_mat(:,ib),p.R.*grids.s.vec + income.meannety1);
            else
                mucnext(:,ib) = aux.utility1(p.risk_aver,coninterp{ib}(p.R*grids.s.vec + income.meannety1))...
                                - p.temptation/(1+p.temptation) * aux.utility1(p.risk_aver,p.R.*grids.s.vec + income.meannety1);
            end
        end
        
        % take expectation over beta
        if numel(p.r) > 1
            emuc = mucnext * heterogeneity.rtrans';
            betastacked = repmat(betagrid',p.nx,p.nb);
        elseif numel(p.risk_aver) > 1
            emuc = mucnext * heterogeneity.ztrans';
            betastacked = repmat(betagrid',p.nx,p.nb);
        else
            emuc = mucnext * heterogeneity.betatrans';
            betastacked = repmat(betagrid',p.nx,1);
        end
        muc1 = (1-p.dieprob) * (1+r_mat) .* betastacked .* emuc ...
                ./ (1+p.savtax*(repmat(grids.s.vec,1,p.nb)>=p.savtaxthresh))...
                + p.dieprob * aux.utility_bequests1(p.bequest_curv,p.bequest_weight,...
                p.bequest_luxury,repmat(grids.s.vec,1,p.nb));
        
        if numel(p.risk_aver) == 1
            con1 = aux.u1inv(p.risk_aver,muc1);
        else
            con1 = aux.u1inv(risk_aver_mat,muc1);
        end
        
        cash1 = con1 + repmat(grids.s.vec,1,p.nb)...
            + p.savtax * max(repmat(grids.s.vec,1,p.nb) - p.savtaxthresh,0);
        
        for ib = 1:p.nb
            savinterp = griddedInterpolant(cash1(:,ib),grids.s.vec,'linear');
            sav(:,ib) = savinterp(grids.x.vec_norisk);
        end
        sav(sav<p.borrow_lim) = p.borrow_lim;
        con = repmat(grids.x.vec_norisk,1,p.nb) - sav...
                                - p.savtax * max(sav - p.savtaxthresh,0);
        
        cdiff = max(abs(con(:)-conlast(:)));
        if mod(iter,100) ==0
            disp([' EGP Iteration (no risk) ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end
    norisk.EGP_cdiff = cdiff;
    norisk.con = con;
    norisk.sav = sav;
    
    % Store consumption and savings function interpolants
    for ib = 1:p.nb
        norisk.coninterp{ib} = griddedInterpolant(grids.x.vec_norisk,con(:,ib),'linear');
        norisk.savinterp{ib} = griddedInterpolant(grids.x.vec_norisk,sav(:,ib),'linear');
    end
    
   
end