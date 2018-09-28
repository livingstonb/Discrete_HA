function [norisk,cdiff] = solve_EGP_deterministic(p,...
    xgrid,sgrid,prefs,income)

    mucnext = zeros(p.nx,p.nb);
    coninterp = cell(1,p.nb);
    sav = zeros(p.nx,p.nb);
    mucnext= zeros(p.nx,p.nb);
    muc1= zeros(p.nx,p.nb);
    con1= zeros(p.nx,p.nb);
    cash1= zeros(p.nx,p.nb);
    
    if  p.nb == 1
        betagrid = p.beta;
    elseif p.nb ==2 
        betagrid = [p.beta-p.betawidth;p.beta+p.betawidth];
    end

    con = p.r * xgrid.norisk_wide;

    iter = 0;
    cdiff = 1000;
    while iter <= p.max_iter && cdiff>p.tol_iter
        iter = iter + 1;
        
        conlast = con;
        muc  = prefs.u1(conlast);
        
        for ib = 1:p.nb
            coninterp{ib} = griddedInterpolant(xgrid.norisk_short,conlast(:,ib),'linear');
            mucnext(:,ib) = prefs.u1(coninterp{ib}(p.R.*sgrid.short + income.meannety));
        end
        
        emuc = mucnext * prefs.betatrans';
        muc1 = (1-p.dieprob) * p.R * repmat(betagrid',p.nx,1) .* emuc ./ (1+p.savtax*(sgrid.norisk_wide>=p.savtaxthresh))...
            + p.dieprob * prefs.beq1(sgrid.norisk_wide);
        
        con1 = prefs.u1inv(muc1);
        cash1 = con1 + sgrid.norisk_wide + p.savtax * max(sgrid.norisk_wide - p.savtaxthresh,0);
        
        sav = interp1(cash1(:),sgrid.norisk_wide,xgrid.norisk_wide(:),'linear','extrap');
        
        % apply borrowing limit
%         idx = bsxfun(@lt,xgrid.norisk_wide,cash1(1,:));
%         sav(idx) = p.borrow_lim;
%         for ib = 1:p.nb
%             idx_ib = idx(:,ib);
%             sav(~idx_ib,ib) = interp1(cash1(~idx_ib,ib),sgrid.short(~idx_ib),xgrid.norisk_wide(~idx_ib,ib),'linear','extrap');
%         end

        sav(sav<p.borrow_lim) = p.borrow_lim;
        

        con = xgrid.norisk_wide - sav - p.savtax * max(sav - p.savtaxthresh,0);
        
        cdiff = max(abs(con(:)-conlast(:)))
        if p.Display >=1 && mod(iter,50) ==0
            disp([' EGP Iteration (no risk) ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end
    
    norisk.con = con;
    norisk.sav = sav;
    
    % Store consumption and savings function interpolants
    for ib = 1:p.nb
        norisk.coninterp{ib} = griddedInterpolant(xgrid.norisk_short,con(:,ib),'linear');
        norisk.savinterp{ib} = griddedInterpolant(xgrid.norisk_short,sav(:,ib),'linear');
    end
    
   
end