function results = egp_AR1_IID_tax_recode(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017

    %% INCOME GRIDS

    %persistent income: rowenhurst
    if p.LoadIncomeProcess == 1
        logyPgrid = load('QuarterlyIncomeDynamics/TransitoryContinuous/logyPgrid.txt');
        yPdist = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPdist.txt');
        yPtrans = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPtrans.txt');
        p.nyP = length(logyPgrid);
        logyPgrid = logyPgrid';
    elseif p.nyP>1
        [logyPgrid, yPtrans, yPdist] = rouwenhorst(p.nyP, -0.5*p.sd_logyP^2, p.sd_logyP, p.rho_logyP);
    else
        logyPgrid = 0;
        yPdist = 1;
        yPtrans = 1;
    end  

    yPgrid = exp(logyPgrid);
    yPcumdist = cumsum(yPdist,1);
    yPcumtrans = cumsum(yPtrans,2);
    
    if size(yPgrid,2)>1
        error('yPgrid is a row vector, must be column vector')
    end
    if size(yPdist,2)>1
        error('yPdist is a row vector, must be column vector')
    end
    

    % transitory income: disretize normal distribution
    if p.LoadIncomeProcess == 1
        p.sig2T = load('QuarterlyIncomeDynamics/TransitoryContinuous/sig2T.txt');
        p.lambdaT = load('QuarterlyIncomeDynamics/TransitoryContinuous/lambdaT.txt');
    end

    if p.nyT>1

        %moments of mixture distribution
        lmu2 = p.lambdaT.*p.sd_logyT^2;
        lmu4 = 3.*p.lambdaT.*(p.sd_logyT^4);

        %fit thjose moments
        optionsNLLS = optimoptions(@lsqnonlin,'Display','Off');
        lpar = lsqnonlin(@(lp)discretize_normal_var_kurt(lp,p.nyT,lmu2,lmu4),[2 0.1],[],[],optionsNLLS);
        [lf,lx,lp] = discretize_normal_var_kurt(lpar,p.nyT,lmu2,lmu4);
        logyTgrid = lx;
        yTdist = lp;
        yTcumdist = cumsum(yTdist,1);

    elseif nyT==1
        logyTgrid = 0;
        yTdist = 1;
    end
    
    yTgrid = exp(logyTgrid);
    
    if size(yTgrid,2)>1
        error('yTgrid is a row vector, must be column vector')
    end
    if size(yTdist,2)>1
        error('yTdist is a row vector, must be column vector')
    end

    % fixed effect
    if p.nyF>1
        width = fzero(@(x)discrete_normal(p.nyF,-0.5*p.sd_logyF^2 ,p.sd_logyF ,x),2);
        [~,logyFgrid,yFdist] = discrete_normal(p.nyF,-0.5*p.sd_logyF^2 ,p.sd_logyF ,width);
    elseif p.nyF==1
        logyFgrid = 0;
        yFdist = 1;
    end
    yFgrid = exp(logyFgrid);
    yFcumdist = cumsum(yFdist,1);

    if size(yFgrid,2)>1
        error('yFgrid is a row vector, must be column vector')
    end
    if size(yFdist,2)>1
        error('yFdist is a row vector, must be column vector')
    end

    % transition probabilities for yP-yF combined grid
    ytrans = kron(eye(p.nyF),yPtrans);

    % length of full xgrid
    p.N = p.nx*p.nyF*p.nyP*p.nb;

    %% DISCOUNT FACTOR

    if p.IterateBeta == 0
        p.maxiterAY = 1;
        % final beta
        beta = p.beta0;
    end

    if p.IterateBeta == 1
        % initial condition for beta iteration
        p.beta0 = (p.betaH + p.betaL)/2;
    end

    %initial discount factor grid
    if  p.nb == 1
        betadist = 1;
        betatrans = 1;
    elseif p.nb ==2 
        betadist = [0.5;0.5];
        betatrans = [1-p.betaswitch p.betaswitch; p.betaswitch 1-p.betaswitch]; %transitions on average once every 40 years;
    else
        error('nb must be 1 or 2');
    end
    betacumdist = cumsum(betadist);
    betacumtrans = cumsum(betatrans,2);


    %% ASSET AND INCOME GRIDS

    sgrid = linspace(0,1,p.nx)';
    sgrid = sgrid.^(1./p.xgrid_par);
    sgrid = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid;

    sgrid_wide = repmat(sgrid,1,p.nyP*p.nyF*p.nb);
    p.ns = p.nx;

    % construct matrix of y combinationsx
    ymat = kron(yPgrid,ones(p.ns,1));
    ymat = repmat(ymat,p.nyF,1) .* kron(yFgrid,ones(p.nyP*p.ns,1));
    ymat = repmat(ymat,p.nb,1)*yTgrid';

    % distribution of ymat (values are repeated nb*nx times)
    ymatdist = kron(yPdist,ones(p.ns,1));
    ymatdist = repmat(ymatdist,p.nyF,1) .* kron(yFdist,ones(p.nyP*p.ns,1));
    ymatdist = repmat(ymatdist,p.nb,1)*yTdist';

    % find mean y
    % isolate unique (yT,yF,yP) combinations
    ymat_yvals = ymat(1:p.ns:p.ns*p.nyF*p.nyP,:);
    ymatdist_pvals = ymatdist(1:p.ns:p.ns*p.nyF*p.nyP,:);
    temp = sortrows([ymat_yvals(:) ymatdist_pvals(:)],1);
    ysortvals = temp(:,1);
    ysortpvals = temp(:,2);
    ycumdist = cumsum(ysortpvals);
    meany = ymat_yvals(:)'*ymatdist_pvals(:);
    original_meany = meany;
    
    % normalize gross income to have mean 1
    if p.NormalizeY == 1
        ymat = ymat/meany;
        ymat_yvals = ymat_yvals/meany;
        ysortvals = ysortvals/meany;
        meany = 1;
    end
    totgrossy = meany;

    % find tax threshold on labor income
    if numel(ysortvals)>1
        labtaxthresh = lininterp1(ycumdist,ysortvals,p.labtaxthreshpc);
    else
        labtaxthresh = 0;
    end    

    % find net income
    totgrossyhigh = max(ymat_yvals(:)-labtaxthresh,0)'*ymatdist_pvals(:);
    lumptransfer = p.labtaxlow*totgrossy + p.labtaxhigh*totgrossyhigh;
    % netymat is N by nyT matrix
    netymat = lumptransfer + (1-p.labtaxlow)*ymat - p.labtaxhigh*max(ymat-labtaxthresh,0);
    netymat_yvals = netymat(1:p.ns:p.ns*p.nyF*p.nyP,:);
    meannety = netymat_yvals(:)'*ymatdist_pvals(:);

    % xgrid, indexed by beta,yF,yP,x (N by 1 matrix)
    % cash on hand grid: different min points for each value of (iyP,iyF)
    xgrid = sgrid_wide(:) + min(netymat,[],2);

    %% UTILITY FUNCTION, BEQUEST FUNCTION

    if p.risk_aver==1
        u = @(c)log(c);
        beq = @(a) p.bequest_weight.* log(a+ p.bequest_luxury);
    else    
        u = @(c)(c.^(1-p.risk_aver)-1)./(1-p.risk_aver);
        beq = @(a) p.bequest_weight.*((a+p.bequest_luxury).^(1-p.risk_aver)-1)./(1-p.risk_aver);
    end    

    u1 = @(c) c.^(-p.risk_aver);
    u1inv = @(u) u.^(-1./p.risk_aver);

    beq1 = @(a) p.bequest_weight.*(a+p.bequest_luxury).^(-p.risk_aver);


    %% Model Solution
    xgrid_wide = reshape(xgrid,p.ns,p.nyP*p.nyF*p.nb);

    if p.IterateBeta == 1
        ergodic_tol = 1e-6;
        iterate_EGP = @(x) solve_EGP(x,p,...
            xgrid_wide,ytrans,betatrans,sgrid_wide,u1,u1inv,netymat,...
            yTdist,beq1,yPtrans,meany,ergodic_tol);

        beta_lb = 1e-3;
        if p.nb == 1
            beta_ub = p.betaH - 1e-5;
        else
            beta_ub = p.betaH - 1e-5 - betawidth;
        end
        % Max fzero iterations set to p.max_evals
        check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.max_evals);
        options = optimset('TolX',1e-6,'OutputFcn',check_evals);
        [beta,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);
        
        if exitflag ~= 1
            results = struct();
            results.issues = {'NoBetaConv'};
            return
        end
        results.beta = beta;
    end

    ergodic_tol = 1e-7;
    [~,con,sav,state_dist,cdiff] = solve_EGP(beta,p,...
        xgrid_wide,ytrans,betatrans,sgrid_wide,u1,u1inv,netymat,...
        yTdist,beq1,yPtrans,meany,ergodic_tol);
    
    % Reshape policy functions for use later
    con_wide = reshape(con,p.nx,p.N/p.nx);
    sav_wide = reshape(sav,p.nx,p.N/p.nx);
    newdim = [p.nx p.nyP p.nyF p.nb];
    con_multidim = reshape(con,newdim);
    sav_multidim = reshape(sav,newdim);
    xgrid_multidim = reshape(xgrid,newdim);
    state_dist_multidim = reshape(state_dist,newdim);
    
    %% Store important moments
    results.mean_s = sav' * state_dist;
    results.mean_x = xgrid' * state_dist;
    results.mean_grossy = (ymat*yTdist)' * state_dist;
    results.mean_loggrossy = (log(ymat)*yTdist)' * state_dist;
    results.mean_nety = (netymat*yTdist)' * state_dist;
    results.mean_lognety = (log(netymat)*yTdist)' * state_dist;
    results.var_loggrossy = state_dist' * (log(ymat) - results.mean_loggrossy).^2 * yTdist;
    results.var_lognety = state_dist' * (log(netymat)- results.mean_lognety).^2 * yTdist;
    
    % Error checks
    mean_x_check = (1+p.r)*results.mean_s + results.mean_nety;
    temp = permute(state_dist_multidim,[2 1 3 4]);
    yPdist_check = sum(sum(sum(temp,4),3),2);
    

    %% Store problems
    results.issues = {};
    if cdiff > p.tol_iter
        results.issues = [results.issues,'NoEGPConv'];
    end
    if (results.mean_s<p.targetAY-1) || (results.mean_s>p.targetAY+1)
        results.issues = [results.issues,'BadAY'];
    end
    if abs((results.mean_x-mean_x_check)/results.mean_x)> 1e-3
        results.issues = [results.issues,'DistNotStationary'];
    end
    if abs((meannety-results.mean_nety)/meannety) > 1e-3
        results.issues = [results.issues,'BadNetIncomeMean'];
    end
    if norm(yPdist_check-yPdist) > 1e-3
        results.issues = [results.issues,'BadGrossIncDist'];
    end
    if min(state_dist) < - 0.01
        results.issues = [results.issues,'NegativeStateProbability'];
    end

    %% WEALTH DISTRIBUTION
    temp = sortrows([sav state_dist]);
    sav_sort = temp(:,1);
    state_dist_sort = temp(:,2);
    state_dist_cum = cumsum(state_dist_sort);
    
    results.frac_constrained = (sav<=p.borrow_lim)' * state_dist;
    results.frac_less5perc_labincome = (sav<0.05)' * state_dist;
    % wealth percentiles;
    percentiles = [0.1 0.25 0.5 0.9 0.99];
    wealthps = zeros(numel(percentiles),1);
    count = 1;
    for percentile = percentiles
        [~,pind] = max(percentile<state_dist_cum,[],1);
        wealthps(count) = sav_sort(pind);
        count = count + 1;
    end
    results.p10wealth = wealthps(1);
    results.p25wealth = wealthps(2);
    results.p50wealth = wealthps(3);
    results.p90wealth = wealthps(4);
    results.p99wealth = wealthps(5);
    
    binwidth = 0.25;
    bins = 0:binwidth:p.xmax;
    values = zeros(p.xmax+1,1);
    ibin = 1;
    for bin = bins
        if bin < p.xmax
            idx = (sav_sort>=bin) & (sav_sort<bin+binwidth);
        else
            idx = (sav_sort>=bin) & (sav_sort<=bin+binwidth);
        end
        values(ibin) = sum(state_dist_sort(idx));
        ibin = ibin + 1;
    end
    
    
    %% Simulate
    if p.Simulate == 1
        [simulations ssim] = simulate(p,yTcumdist,yFcumdist,...
    yPcumdist,yPcumtrans,yPgrid,yFgrid,yTgrid,labtaxthresh,con_multidim,sav_multidim,xgrid_multidim,...
    lumptransfer,betacumdist,betacumtrans,original_meany);
    else
        simulations =[];
    end

    %% MAKE PLOTS

    if p.MakePlots ==1 

     figure(1);

        %plot for median fixed effect
        if mod(p.nyF,2)==1
            iyF = (p.nyF+1)/2;
        else
            iyF = p.nyF/2;
        end

        % plot for first beta
        iyb = 1;
        % if nb = 1, force plot of first beta
        if p.nb == 1
            iyb = 1;
        end

        % consumption policy function
        subplot(2,4,1);
        plot(xgrid_multidim(:,1,iyF,iyb),con_multidim(:,1,iyF,iyb),'b-',xgrid_multidim(:,p.nyP,iyF,iyb),con_multidim(:,p.nyP,iyF,iyb),'r-','LineWidth',1);
        grid;
        xlim([p.borrow_lim p.xmax]);
        title('Consumption Policy Function');
        legend('Lowest income state','Highest income state');

        % savings policy function
        subplot(2,4,2);
        plot(xgrid_multidim(:,1,iyF),sav_multidim(:,1,iyF)./xgrid_multidim(:,1,iyF),'b-',xgrid_multidim(:,p.nyP,iyF),sav_multidim(:,p.nyP,iyF)./xgrid_multidim(:,p.nyP,iyF),'r-','LineWidth',1);
        hold on;
        plot(sgrid,ones(p.nx,1),'k','LineWidth',0.5);
        hold off;
        grid;
        xlim([p.borrow_lim p.xmax]);
        title('Savings Policy Function s/x');

        % consumption policy function: zoomed in
        subplot(2,4,3);
        plot(xgrid_multidim(:,1,iyF),con_multidim(:,1,iyF),'b-',xgrid_multidim(:,p.nyP,iyF),con_multidim(:,p.nyP,iyF),'r-','LineWidth',2);
        grid;
        xlim([0 4]);
        title('Consumption: Zoomed');

         % savings policy function: zoomed in
        subplot(2,4,4);
        plot(xgrid_multidim(:,1,iyF),sav_multidim(:,1,iyF)./xgrid_multidim(:,1,iyF),'b-',xgrid_multidim(:,p.nyP,iyF),sav_multidim(:,p.nyP,iyF)./xgrid_multidim(:,p.nyP,iyF),'r-','LineWidth',2);
        hold on;
        plot(sgrid,ones(p.nx,1),'k','LineWidth',0.5);
        hold off;
        grid;
        xlim([0 4]);
        title('Savings (s/x): Zoomed');
        
         % gross income distribution
        subplot(2,4,5);
        b = bar(ysortvals,ysortpvals);
        b.FaceColor = 'blue';
        b.EdgeColor = 'blue';
        grid;
        xlim([0 10]);
        title('Gross Income PMF');
        
         % asset distribution
        subplot(2,4,6);
        b = bar(bins,values);
        b.FaceColor = 'blue';
        b.EdgeColor = 'blue';
        grid;
        xlim([-0.4 10]);
        ylim([0 1]);
        title('Asset PMF, Binned');

         % simulation convergence
        if p.Simulate == 1
            subplot(2,4,7);
            plot(1:p.Tsim,mean(ssim),'b','LineWidth',2);
            grid;
            xlim([0 p.Tsim]);
            title('Mean savings (sim)');
        end
    end

    %% COMPUTE MPCs
    if p.ComputeMPC ==1
        %theoretical mpc lower bound
        mpclim = p.R*((beta*p.R)^-(1./p.risk_aver))-1;
        Nmpcamount = numel(p.mpcfrac);
        %mpc amounts
        for im = 1:Nmpcamount
            mpcamount{im} = p.mpcfrac{im} * meany;
            xgrid_mpc_wide{im} = xgrid_wide + mpcamount{im};
        end

        results.mpcamount = mpcamount;
        
        % Create interpolants for computing mpc's
        for col = 1:p.N/p.nx
                coninterp{col} = griddedInterpolant(xgrid_wide(:,col),con_wide(:,col),'linear');
        end
        
        % mpc functions
        for im = 1:Nmpcamount
            mpc{im} = zeros(p.nx,p.N/p.nx);
            % iterate over (yP,yF,beta)
            for col = 1:p.N/p.nx 
                mpc{im}(:,col) = (coninterp{col}(xgrid_mpc_wide{im}(:,col))...
                    - con_wide(:,col))/mpcamount{im};     
            end
            % average mpc
            results.avg_mpc{im} = mpc{im}(:)' * state_dist;
        end
    end
    
%% Print Results
    if p.PrintStats == 1
        print_statistics(results,simulations,p);
    end
end