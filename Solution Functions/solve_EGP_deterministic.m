function norisk = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income,direct_results)
    % This function uses the method of endogenous grid points to find the
    % policy functions of the deterministic model. Output is in the
    % 'norisk' structure.

    mucnext = zeros(p.nx,p.nb);
    coninterp = cell(1,p.nb);
    sav = zeros(p.nx,p.nb);
    mucnext= zeros(p.nx,p.nb);
    muc1= zeros(p.nx,p.nb);
    con1= zeros(p.nx,p.nb);
    cash1= zeros(p.nx,p.nb);
    
    betagrid = direct_results.beta + prefs.betagrid0;
    
    if numel(p.risk_aver) > 1
        risk_aver_mat = repmat(p.risk_aver,p.nx,1);
    end

    con = p.r * repmat(sgrid.short,1,p.nb) + income.meany1;

    iter = 0;
    cdiff = 1000;
    while iter <= p.max_iter && cdiff>p.tol_iter
        iter = iter + 1;
        
        conlast = con;
        
        for ib = 1:p.nb
            coninterp{ib} = griddedInterpolant(xgrid.norisk_short,conlast(:,ib),'linear');
            
            % cash-on-hand is just Rs + meany
            if numel(p.risk_aver) == 1
                mucnext(:,ib) = prefs.u1(coninterp{ib}(p.R.*sgrid.short + income.meany1))...
                                - p.temptation/(1+p.temptation) * prefs.u1(p.R.*sgrid.short + income.meany1);
            else
                mucnext(:,ib) = prefs.u1(risk_aver_mat(:,ib),coninterp{ib}(p.R.*sgrid.short + income.meany1))...
                                - p.temptation/(1+p.temptation) * prefs.u1(risk_aver_mat(:,ib),p.R.*sgrid.short + income.meany1);
            end
        end
        
        % take expectation over beta
        if numel(p.risk_aver) == 1
            emuc = mucnext * prefs.betatrans';
            betastacked = repmat(betagrid',p.nx,1);
        else
            emuc = mucnext * prefs.IEStrans';
            betastacked = repmat(betagrid',p.nx,p.nb);
        end
        muc1 = (1-p.dieprob) * p.R * betastacked .* emuc ...
                ./ (1+p.savtax*(repmat(sgrid.short,1,p.nb)>=p.savtaxthresh))...
                + p.dieprob * prefs.beq1(repmat(sgrid.short,1,p.nb));
        
        if numel(p.risk_aver) == 1
            con1 = prefs.u1inv(muc1);
        else
            con1 = prefs.u1inv(risk_aver_mat,muc1);
        end
        
        cash1 = con1 + repmat(sgrid.short,1,p.nb)...
            + p.savtax * max(repmat(sgrid.short,1,p.nb) - p.savtaxthresh,0);
        
        for ib = 1:p.nb
            savinterp = griddedInterpolant(cash1(:,ib),sgrid.short,'linear');
            sav(:,ib) = savinterp(xgrid.norisk_short);
        end
        sav(sav<p.borrow_lim) = p.borrow_lim;
        con = repmat(xgrid.norisk_short,1,p.nb) - sav...
                                - p.savtax * max(sav - p.savtaxthresh,0);
        
        cdiff = max(abs(con(:)-conlast(:)));
        if p.Display ==1 && mod(iter,100) ==0
            disp([' EGP Iteration (no risk) ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end
    norisk.EGP_cdiff = cdiff;
    norisk.con = con;
    norisk.sav = sav;
    
    % Store consumption and savings function interpolants
    for ib = 1:p.nb
        norisk.coninterp{ib} = griddedInterpolant(xgrid.norisk_short,con(:,ib),'linear');
        norisk.savinterp{ib} = griddedInterpolant(xgrid.norisk_short,sav(:,ib),'linear');
    end
    
   
end