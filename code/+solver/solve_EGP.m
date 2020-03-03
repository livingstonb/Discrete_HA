function model = solve_EGP(p, grids, heterogeneity,...
    income, nextmpcshock, prevmodel)
    % This function performs the method of endogenous grid points to find
    % saving and consumption policy functions. It also calls 
    % find_stationary() to compute the stationary distribution over states 
    % via direct methods (rather than simulations) and stores the results 
    % in the 'model' structure.
    %
    % To compute the MPCs out of news, it is necessary for the policy function
    % to reflect the expectation of a future shock. For these cases,
    % the policy functions in 'prevmodel' are used. The variable 'nextmpcshock'
    % is nonzero when a shock is expected next period.
    %
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    %% ----------------------------------------------------
    % USEFUL OBJECTS/ARRAYS
    % -----------------------------------------------------
    ss_dims = [p.nx, p.nyP, p.nyF, p.nb];
    ss_dims_aug = [ss_dims p.nyT];

    repmat_to_state_space = ...
        @(arr) aux.repmat_auto(arr, ss_dims);
    repmat_to_state_space_aug = ...
        @(arr) aux.repmat_auto(arr, ss_dims_aug);

    reshape_to_state_space = ...
        @(arr) reshape(arr, ss_dims);

    sgrid_tax = p.compute_savtax(grids.s.vec);

    r_bc = heterogeneity.r_broadcast;
    R_bc = heterogeneity.R_broadcast;

    tempt_bc = heterogeneity.temptation_broadcast;
    tempt_expr = tempt_bc ./ (1 + tempt_bc);

    betagrid_bc = heterogeneity.betagrid_broadcast;

    % Find xprime as a function of s
    tmp = R_bc .* grids.s.matrix + income.netymatEGP + nextmpcshock;
    xprime_s = repmat_to_state_space_aug(tmp);

    %% ----------------------------------------------------
    % REGION WHERE NEXT PERIOD'S ASSETS GUARANTEED NON-NEG
    % ----------------------------------------------------- 
    % min_nety = min(income.netymat(:));
%     svalid = heterogeneity.R_broadcast .* grids.s.matrix...
%         + min_nety + nextmpcshock >= grids.x.matrix(1,:,:,:);
    
    %% ----------------------------------------------------
    % CONSTRUCT EXPECTATIONS MATRIX, ETC...
    % -----------------------------------------------------
    
    Emat = kron(income.ytrans_live, speye(p.nx));

    % Initial guess for consumption function
    tempt_adj = 0.5 * (max(p.temptation) > 0.05);
    r_mat_adj = max(r_bc, 0.001);
    con = (r_mat_adj + tempt_adj) .* grids.x.matrix;
    con = con(:);
    con(con<=0) = min(con(con>0));
    con = reshape_to_state_space(con);

    %% ----------------------------------------------------
    % EGP ITERATION
    % ----------------------------------------------------- 
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

        % c(x')
        c_xp = get_c_xprime(p, grids, xprime_s, prevmodel, conlast, nextmpcshock);

        % MUC in current period, from Euler equation
        muc_s = get_marginal_util_cons(...
            p, income, grids, c_xp, xprime_s, R_bc,...
            Emat, betagrid_bc, heterogeneity.risk_aver_broadcast,...
            tempt_expr);
     
        % c(s)
        con_s = aux.u1inv(heterogeneity.risk_aver_broadcast, muc_s);
        
        % x(s) = s + stax + c(s)
        x_s = grids.s.vec + sgrid_tax + con_s;

        % interpolate from x(s) to get s(x)
        sav = get_saving_policy(p, grids, x_s);
        sav_tax = p.compute_savtax(sav);

        % updated consumption function, column vec length of
        % length p.nx*p.nyP*p.nyF*p.nb
        conupdate = grids.x.matrix - sav - sav_tax;

        cdiff = max(abs(conupdate(:)-conlast(:)));
        if mod(iter,50) ==0
            disp(['  EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end
    
    if cdiff>p.tol_iter
        % EGP did not converge, don't find stationary distribution
        AYdiff = 100000;
        model.EGP_cdiff = cdiff;
        return
    end

    model.sav = sav;
    model.con = conupdate;
    model.EGP_cdiff = cdiff;

    % create interpolants from optimal policy functions
    % and find saving values associated with xvals
    model.savinterp = cell(p.nyP,p.nyF,p.nb);
    model.coninterp = cell(p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        model.savinterp{iyP,iyF,ib} =  griddedInterpolant(...
            grids.x.matrix(:,iyP,iyF,ib), model.sav(:,iyP,iyF,ib), 'linear');
        model.coninterp{iyP,iyF,ib} = griddedInterpolant(...
            grids.x.matrix(:,iyP,iyF,ib), model.con(:,iyP,iyF,ib), 'linear');

        xmin = grids.x.matrix(1,iyP,iyF,ib);
        model.coninterp_xprime{iyP,iyF,ib} = @(xprime) new_con_interp(...
            model.coninterp{iyP,iyF,ib}, xprime, nextmpcshock, xmin);
    end
    end
    end
end

function result = new_con_interp(old_con_interp, xprime, nextmpcshock, xmin)
    valid = xprime + nextmpcshock >= xmin;
    result = zeros(size(xprime));
    result(valid) = old_con_interp(xprime(valid));
    result(~valid) = 1e-8;
end

function c_xprime = get_c_xprime(p, grids, xp_s, prevmodel, conlast, nextmpcshock)
	% find c as a function of x'
	c_xprime = zeros(size(xp_s));
    dims_nx_nyT = [p.nx, 1, 1, 1, p.nyT];
    
    % if nextmpcshock >= 0
    %     ixp_valid = 1:p.nx;
    % else
    %     xpvalid = grids.x.matrix + nextmpcshock >= grids.x.matrix(1,:,:,:);
    % end

    ixp_valid = 1:p.nx;

	for ib  = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        % if nextmpcshock < 0
        %     ixp_valid = find(xpvalid(:,iyP,iyF,ib));
        % end
    	xp_s_ib_iyF_iyP = xp_s(ixp_valid,iyP,iyF,ib,:);

        if isempty(prevmodel)
            % usual method of EGP
            coninterp = griddedInterpolant(...
                grids.x.matrix(ixp_valid,iyP,iyF,ib), conlast(ixp_valid,iyP,iyF,ib), 'linear');
            c_xprime(ixp_valid,iyP,iyF,ib,:) = reshape(...
                coninterp(xp_s_ib_iyF_iyP(:)), dims_nx_nyT);
        else
            % need to compute IMPC(s,t) for s > 1, where IMPC(s,t) is MPC in period t out of period
            % s shock that was learned about in period 1 < s
            c_xprime(ixp_valid,iyP,iyF,ib,:) = reshape(...
                prevmodel.coninterp_xprime{iyP,iyF,ib}(xp_s_ib_iyF_iyP(:)), dims_nx_nyT);
        end
    end
    end
    end
end

function muc_s = get_marginal_util_cons(...
	p, income, grids, c_xp, xp_s, R_bc,...
    Emat, betagrid_bc, risk_aver_bc, tempt_expr)

    savtaxrate = (1 + p.savtax .* (grids.s.matrix >= p.savtaxthresh));
    savtaxrate = repmat(savtaxrate, [1, 1, 1, p.nb]);

	% First get marginal utility of consumption next period
	muc_c = aux.utility1(risk_aver_bc, c_xp);
    muc_tempt = -tempt_expr .* aux.utility1(risk_aver_bc, xp_s);
    mucnext = reshape(muc_c(:) + muc_tempt(:), [], p.nyT);

    % Integrate
    expectation = Emat * mucnext * income.yTdist;
    expectation = reshape(expectation, [p.nx, p.nyP, p.nyF, p.nb]);

    muc_c_today = R_bc .* betagrid_bc .* expectation;
    muc_beq = aux.utility_bequests1(p.bequest_curv, p.bequest_weight,...
        p.bequest_luxury, grids.s.matrix);

    muc_s = (1 - p.dieprob) * muc_c_today ./ savtaxrate ...
        + p.dieprob * muc_beq;
end

function sav = get_saving_policy(p, grids, x_s)
	% finds s(x), the saving policy function on the
	% cash-on-hand grid

	sav = zeros(size(x_s));
    for ib  = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        for i = 1:p.nx
            if x_s(i,iyP,iyF,ib) < x_s(i+1,iyP,iyF,ib)
                is_valid = i:p.nx;
                break;
            end
        end
        savinterp = griddedInterpolant(x_s(is_valid,iyP,iyF,ib),...
            grids.s.matrix(is_valid,iyP,iyF), 'linear');
        sav(:,iyP,iyF,ib) = savinterp(grids.x.matrix(:,iyP,iyF,ib));
        
        adj = grids.x.matrix(:,iyP,iyF,ib) < x_s(1,iyP,iyF,ib);
        sav(adj,iyP,iyF,ib) = p.borrow_lim;
    end
    end
    end

    % % deal with borrowing limit
    % sav(sav<p.borrow_lim) = p.borrow_lim;
end