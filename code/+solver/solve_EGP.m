function model = solve_EGP(beta,p,grids,heterogeneity,...
    income,nextmpcshock,prevmodel)
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
    % REGION WHERE NEXT PERIOD'S ASSETS GUARANTEED NON-NEG
    % ----------------------------------------------------- 
    min_nety = min(income.netymat(:));
    xvalid = p.R * grids.x.matrix + min_nety + nextmpcshock >= 0;
    svalid = p.R * grids.s.matrix + min_nety + nextmpcshock >= 0;
    
    %% ----------------------------------------------------
    % CONSTRUCT EXPECTATIONS MATRIX, ETC...
    % -----------------------------------------------------
    if numel(p.beta_grid_forced) == 0
        betagrid = beta + heterogeneity.betagrid0;
    else
        betagrid = p.beta_grid_forced;
    end
    
    if p.IterateBeta == 1
        msg = sprintf(' %3.3f',betagrid);
        disp([' Trying betagrid =' msg])
    end

    % Expectations operator (conditional on yT)
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if numel(p.r) > 1
        Emat = kron(heterogeneity.rtrans,kron(income.ytrans,speye(p.nx)));
        r_col = kron(p.r',ones(p.nx*p.nyP*p.nyF,1));
        r_mat = reshape(r_col,[p.nx,p.nyP,p.nyF,numel(p.r)]);
    elseif numel(p.risk_aver) > 1
        Emat = kron(heterogeneity.ztrans,kron(income.ytrans,speye(p.nx)));
        risk_aver_col = kron(p.risk_aver',ones(p.nx*p.nyP*p.nyF,1));
        r_mat = p.r;
    else
        Emat = kron(heterogeneity.betatrans,kron(income.ytrans,speye(p.nx)));
        r_mat = p.r;
    end

    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    if p.temptation > 0.05
        extra = 0.5;
    else
        extra = 0;
    end
    
    con = (r_mat(:) .* (r_mat(:)>=0.001) + 0.001 * (r_mat(:)<0.001) + extra) ...
    	.* repmat(grids.x.matrix(:), p.nb, 1);
    con(con<=0) = min(con(con>0));

    % discount factor matrix, 
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    if (numel(p.risk_aver) > 1) || (numel(p.r) > 1)
        % IES heterogeneity or returns heterogeneity - nb is number of IES or r values
        % betagrid is just beta
        betastacked = speye(p.nyP*p.nyF*p.nx*p.nb) * betagrid;
    else
        % beta heterogeneity
        betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));
        betastacked = sparse(diag(betastacked));
    end

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
        
        % c(x)
        conlast = reshape(conlast,[p.nx p.nyP p.nyF p.nb]);
        
        % x'(s)
        xp_s = get_xprime_s(p,income,grids,r_mat,nextmpcshock);

        % c(x')
        c_xp = get_c_xprime(p,grids,xp_s,prevmodel,conlast,nextmpcshock);
        
        % reshape to take expecation over yT first
        c_xp = reshape(c_xp,[],p.nyT);
        xp_s = reshape(xp_s,[],p.nyT);

        % MUC in current period, from Euler equation
        muc_s = get_marginal_util_cons(...
        	p,income,grids,c_xp,xp_s,r_mat,Emat,betastacked);
     
        % c(s)
        if numel(p.risk_aver) == 1
            con_s = aux.u1inv(p.risk_aver,muc_s);
        else
            con_s = aux.u1inv(risk_aver_col,muc_s);
        end
        
        % x(s) = s + stax + c(s)
        x_s = repmat(grids.s.matrix(:),p.nb,1)...
                        + p.savtax * max(repmat(grids.s.matrix(:),p.nb,1)-p.savtaxthresh,0)...
                        + con_s;
        x_s = reshape(x_s,[p.nx p.nyP p.nyF p.nb]);

        % interpolate from x(s) to get s(x)
        sav = get_saving_policy(p,grids,x_s,svalid,nextmpcshock);

        % updated consumption function, column vec length of
        % length p.nx*p.nyP*p.nyF*p.nb
        conupdate = repmat(grids.x.matrix(:),p.nb,1) - sav(:)...
                            - p.savtax * max(sav(:)-p.savtaxthresh,0);

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
        model.coninterp_xprime{iyP,iyF,ib} = @(xprime) new_con_interp(model.coninterp{iyP,iyF,ib},xprime,nextmpcshock);
    end
    end
    end
end

function result = new_con_interp(old_con_interp, xprime, nextmpcshock)
    valid = xprime + nextmpcshock >= 0;
    result = zeros(size(xprime));
    result(valid) = old_con_interp(xprime(valid));
    result(~valid) = 1e-8;
end
    

function xprime_s = get_xprime_s(p,income,grids,r_mat,nextmpcshock)
	% find xprime as a function of s

    xprime_s = (1+r_mat) .* grids.s.matrix + income.netymatEGP + nextmpcshock;
end

function c_xprime = get_c_xprime(p,grids,xp_s,prevmodel,conlast,nextmpcshock)
	% find c as a function of x'
	c_xprime = zeros(p.nx,p.nyP,p.nyF,p.nb,p.nyT);
    
    if nextmpcshock >= 0
        ixp_valid = 1:p.nx;
    else
        xpvalid = grids.x.matrix + nextmpcshock >= 0;
    end

	for ib  = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        if nextmpcshock < 0
            ixp_valid = find(xpvalid(:,iyP,iyF));
        end
    	xp_s_ib_iyF_iyP = xp_s(ixp_valid,iyP,iyF,ib,:);
        if isempty(prevmodel)
            % usual method of EGP
            coninterp = griddedInterpolant(grids.x.matrix(ixp_valid,iyP,iyF),conlast(ixp_valid,iyP,iyF,ib),'linear');
            c_xprime(ixp_valid,iyP,iyF,ib,:) = reshape(coninterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
        else
            % need to compute IMPC(s,t) for s > 1, where IMPC(s,t) is MPC in period t out of period
            % s shock that was learned about in period 1 < s
            c_xprime(ixp_valid,iyP,iyF,ib,:) = reshape(prevmodel.coninterp_xprime{iyP,iyF,ib}(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
        end
    end
    end
    end
end

function muc_s = get_marginal_util_cons(...
	p,income,grids,c_xp,xp_s,r_mat,Emat,betastacked)

	% first get marginal utility of consumption next period
	if numel(p.risk_aver) > 1
		risk_aver_col = kron(p.risk_aver',ones(p.nx*p.nyP*p.nyF,1));
        risk_aver_col_yT = repmat(risk_aver_col,1,p.nyT);
        mucnext = aux.utility1(risk_aver_col_yT,c_xp)...
            - p.temptation/(1+p.temptation) * aux.utility1(risk_aver_col_yT,xp_s);
    else
        mucnext = aux.utility1(p.risk_aver,c_xp) ...
            - p.temptation/(1+p.temptation) * aux.utility1(p.risk_aver,xp_s);
    end

    % now get MUC this period as a function of s
    savtaxrate  = (1+p.savtax.*(repmat(grids.s.matrix(:),p.nb,1)>=p.savtaxthresh));
    mu_consumption = (1+r_mat(:)).*betastacked*Emat*(mucnext*income.yTdist);
    mu_bequest = aux.utility_bequests1(p.bequest_curv,p.bequest_weight,...
                    p.bequest_luxury,repmat(grids.s.matrix(:),p.nb,1));
    muc_s = (1-p.dieprob) * mu_consumption ./ savtaxrate...
                                            + p.dieprob * mu_bequest;
end

function sav = get_saving_policy(p, grids, x_s, svalid, nextmpcshock)
	% finds s(x), the saving policy function on the
	% cash-on-hand grid

	sav = zeros(p.nx,p.nyP,p.nyF,p.nb);
    for ib  = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        for i = 1:p.nx
            if (x_s(i,iyP,iyF,ib) < x_s(i+1,iyP,iyF,ib)) && svalid(i,1,1,1)
                is_valid = i:p.nx;
                break;
            end
        end
        savinterp = griddedInterpolant(x_s(is_valid,iyP,iyF,ib),grids.s.matrix(is_valid,iyP,iyF),'linear');
        sav(:,iyP,iyF,ib) = savinterp(grids.x.matrix(:,iyP,iyF)); 
    end
    end
    end

    % deal with borrowing limit
    sav(sav<p.borrow_lim) = p.borrow_lim;
end