function norisk = solve_EGP_deterministic(p, grids,...
    heterogeneity, income, direct_results)
    % This function uses the method of endogenous grid points to find the
    % policy functions of the deterministic model. Output is in the
    % 'norisk' structure.
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    coninterp = cell(1,p.nb);
    sav = zeros(p.nx,p.nb);
    mucnext= zeros(p.nx,p.nb);

    sgrid_repeated = repmat(grids.s.vec, p.nb, 1);
    sgrid_tax = p.compute_savtax(sgrid_repeated);

    % if numel(p.r) > 1
    %     Emat = kron(heterogeneity.rtrans, kron(income.ytrans, speye(p.nx)));
    %     r_col = kron(p.r', ones(p.nx, 1));
    %     r_mat = reshape(r_col, [p.nx, numel(p.r)]);
    % else
    %     r_mat = p.r;
    % end

    if numel(p.r) > 1
        exog_trans = heterogeneity.rtrans;
    elseif numel(p.nbeta) > 1
        exog_trans = heterogeneity.betatrans;
    else
        exog_trans = heterogeneity.ztrans;
    end

    Emat = kron(exog_trans, speye(p.nx));
    r_mat = reshape(p.r, [1 numel(p.r)]);

    if numel(p.risk_aver) > 1
        risk_aver_mat = kron(p.risk_aver, ones(p.nx,1));
    end

    % initial guess for consumption function, stacked state combinations
    if p.temptation > 0.005
        extra = 0.5;
    else
        extra = 0;
    end
    
    con = (r_mat .* (r_mat>=0.001) + 0.001 * (r_mat<0.001) + extra) .* grids.x.matrix_norisk;
    con_vec = con(:);
    con(con<=0) = min(con_vec(con_vec>0));

    iter = 0;
    cdiff = 1000;
    while iter <= p.max_iter && cdiff>p.tol_iter
        iter = iter + 1;
        
        conlast = con;
        
        for ib = 1:p.nb

            coninterp{ib} = griddedInterpolant(grids.x.matrix_norisk(:,ib), conlast(:,ib), 'linear');
            if numel(p.r) > 1
                con_ib = coninterp{ib}(p.R(ir) * grids.s.vec + income.meannety1);
            else
                con_ib = coninterp{ib}(p.R * grids.s.vec + income.meannety1);
            end
            
            % cash-on-hand is just Rs + meany
            if numel(p.r) > 1
                mucnext(:,ib) = aux.utility1(p.risk_aver, con_ib)...
                                - p.temptation/(1+p.temptation) * aux.utility1(p.risk_aver, con_ib);
            elseif numel(p.risk_aver) > 1
                mucnext(:,ib) = aux.utility1(risk_aver_mat(:,ib), con_ib)...
                                - p.temptation/(1+p.temptation) * aux.utility1(risk_aver_mat(:,ib), con_ib);
            elseif numel(p.temptation) > 1
                mucnext(:,ib) = aux.utility1(p.risk_aver, con_ib)...
                                - p.temptation(ib) / (1+p.temptation(ib)) * aux.utility1(p.risk_aver, con_ib);
            else
                mucnext(:,ib) = aux.utility1(p.risk_aver, con_ib)...
                                - p.temptation/(1+p.temptation) * aux.utility1(p.risk_aver, con_ib);
            end
        end
        
        % take expectation over beta
        betastacked = kron(heterogeneity.betagrid(:), ones(p.nx,1));
        if numel(p.nbeta) == 1
            betastacked = repmat(betastacked, p.nb, 1);
        end

        % if numel(p.r) > 1
        %     emuc = mucnext * transpose(heterogeneity.rtrans);
        % elseif numel(p.nbeta) > 1
        %     emuc = mucnext * transpose(heterogeneity.betatrans);
        % else
        %     emuc = mucnext * transpose(heterogeneity.ztrans);
        % end

        emuc = Emat * mucnext(:);

        muc1 = (1-p.dieprob) * (1+r_mat) .* betastacked .* emuc ...
                ./ (1+p.savtax*(repmat(grids.s.vec,1,p.nb)>=p.savtaxthresh))...
                + p.dieprob * aux.utility_bequests1(p.bequest_curv,p.bequest_weight,...
                p.bequest_luxury,repmat(grids.s.vec,1,p.nb));
        
        if numel(p.risk_aver) == 1
            con1 = aux.u1inv(p.risk_aver, muc1);
        else
            con1 = aux.u1inv(risk_aver_mat, muc1);
        end
        
        cash1 = con1 + repmat(grids.s.vec, 1, p.nb) + sgrid_tax;
        
        for ib = 1:p.nb
            savinterp = griddedInterpolant(cash1(:,ib), grids.s.vec, 'linear');
            sav(:,ib) = max(savinterp(grids.x.matrix_norisk(:,ib)), p.borrow_lim);

            con = grids.x.matrix_norisk(:,ib) - sav(:,ib) - p.compute_savtax(sav(:,ib));
        end

        
        
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
        norisk.coninterp{ib} = griddedInterpolant(grids.x.matrix_norisk(:,ib), con(:,ib), 'linear');
        norisk.savinterp{ib} = griddedInterpolant(grids.x.matrix_norisk(:,ib), sav(:,ib), 'linear');
    end
    
   
end