function norisk = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income)
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
    
    betagrid = p.beta + prefs.betagrid0;

    con = p.r * repmat(sgrid.short,1,p.nb) + income.meannety;

    iter = 0;
    cdiff = 1000;
    while iter <= p.max_iter && cdiff>p.tol_iter
        iter = iter + 1;
        
        conlast = con;
        
        for ib = 1:p.nb
            coninterp{ib} = griddedInterpolant(xgrid.norisk_short,conlast(:,ib),'linear');
            % cash-on-hand is just Rs + meannety
            mucnext(:,ib) = prefs.u1(coninterp{ib}(p.R.*sgrid.short + income.meannety));
        end
        
        % take expectation over beta
        emuc = mucnext * prefs.betatrans';
        muc1 = (1-p.dieprob) * p.R * repmat(betagrid',p.nx,1) .* emuc ...
                ./ (1+p.savtax*(repmat(sgrid.short,1,p.nb)>=p.savtaxthresh))...
                + p.dieprob * prefs.beq1(repmat(sgrid.short,1,p.nb));
        
        con1 = prefs.u1inv(muc1);
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