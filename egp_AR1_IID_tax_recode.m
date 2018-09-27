function results = egp_AR1_IID_tax_recode(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    %% ADJUST PARAMETERS FOR DATA FREQUENCY, OTHER CHOICES

    p.r = p.r/p.freq;
    p.R = 1 + p.r;
    p.beta0 = p.beta0^(1/p.freq);
    p.dieprob = p.dieprob^(1/p.freq);
    p.betaswitch = p.betaswitch^(1/p.freq);
    p.beta = p.betaL^(1/p.freq);
    p.betaH = 1/((p.R)*(1-p.dieprob));
    

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

    sgrid.orig = linspace(0,1,p.nx)';
    sgrid.orig = sgrid.orig.^(1./p.xgrid_par);
    sgrid.orig = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid.orig;
    sgrid.short = sgrid.orig;

    sgrid.wide = repmat(sgrid.short,[1 p.nyP p.nyF p.nb]);
    p.ns = p.nx;

    % construct matrix of y combinationsx
    ymat = repmat(yPgrid,p.nyF,1) .* kron(yFgrid,ones(p.nyP,1)) * yTgrid';

    % distribution of ymat (values are repeated nb*nx times)
    ymatdist = repmat(yPdist,p.nyF,1) .* kron(yFdist,ones(p.nyP,1)) * yTdist';

    % find mean y
    % isolate unique (yT,yF,yP) combinations
    temp = sortrows([ymat(:) ymatdist(:)],1);
    ysort = temp(:,1);
    ysortdist = temp(:,2);
    ycumdist = cumsum(ysortdist);
    meany = ymat(:)'*ymatdist(:);
    original_meany = meany;
    
    % normalize gross income to have ap.nxlongual mean 1
    if p.NormalizeY == 1
        ymat = ymat/(meany*p.freq);
        ysort = ysort/(meany*p.freq);
        meany = 1/p.freq;
    end
    totgrossy = meany;

    % find tax threshold on labor income
    if numel(ysort)>1
        labtaxthresh = lininterp1(ycumdist,ysort,p.labtaxthreshpc);
    else
        labtaxthresh = 0;
    end    

    % find net income
    totgrossyhigh = max(ymat(:)-labtaxthresh,0)'*ymatdist(:);
    lumptransfer = p.labtaxlow*totgrossy + p.labtaxhigh*totgrossyhigh;
    % netymat is N by nyT matrix
    netymat = lumptransfer + (1-p.labtaxlow)*ymat - p.labtaxhigh*max(ymat-labtaxthresh,0);
    meap.nxlongety = netymat(:)'*ymatdist(:);

    % xgrid, indexed by beta,yF,yP,x (N by 1 matrix)
    % cash on hand grid: different min points for each value of (iyP,iyF)
    xgrid.orig = sgrid.wide(:) + min(kron(netymat,ones(p.nx,1)),[],2);
    
    % Store income variables in a structure
    newfields = {'ymat','netymat','meany','original_meany','yPgrid',...
        'yTgrid','yFgrid','yPdist','yTdist','yFdist','yPcumtrans',...
        'yPtrans','yPcumdist','yFcumdist','yTcumdist','ytrans'};
    for i = 1:numel(newfields)
        income.(newfields{i}) = eval(newfields{i});
    end

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

    %% ORGANIZE VARIABLES
    % Reshape policy functions for use later
    xgrid.orig_wide = reshape(xgrid.orig,[p.nx p.nyP p.nyF p.nb]);

    %% MODEL SOLUTION

    if p.IterateBeta == 1
        if p.FastIter == 1
            [beta,exitflag] = iterate_beta(p,...
            xgrid,sgrid,betatrans,u1,beq1,u1inv,income);
        else
            ergodic_tol = 1e-5;
            gridsize = 100;
            iterate_EGP = @(x) solve_EGP(x,p,...
                xgrid,sgrid,betatrans,u1,beq1,u1inv,ergodic_tol,income,ExpandGrid);

            beta_lb = p.betaL;
            if p.nb == 1
                beta_ub = p.betaH - 1e-5;
            else
                beta_ub = p.betaH - 1e-5 - betawidth;
            end

            check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
            options = optimset('TolX',p.tolAY,'OutputFcn',check_evals);
            [beta,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);
        end

    if exitflag ~= 1
        results = struct();
        results.issues = {'NoBetaConv'};
        return
    end
    results.beta = beta;
    end

    ergodic_tol = 1e-6;
    [~,con.orig,sav.orig,con.final,sav.final,state_dist.final,cdiff,xgrid.final] = solve_EGP(beta,p,...
        xgrid,sgrid,betatrans,u1,beq1,u1inv,ergodic_tol,income,p.nxlong);
    
    %% Store important moments
    
    ymat_onxgrid = kron(ymat,ones(p.nxlong,1));
    netymat_onxgrid = kron(netymat,ones(p.nxlong,1));
    
    results.mean_s = sav.final' * state_dist.final;
    results.mean_x = xgrid.final(:)' * state_dist.final;
    results.mean_grossy = (ymat_onxgrid*yTdist)' * state_dist.final;
    results.mean_loggrossy = (log(ymat_onxgrid)*yTdist)' * state_dist.final;
    results.mean_nety = (netymat_onxgrid*yTdist)' * state_dist.final;
    results.mean_lognety = (log(netymat_onxgrid)*yTdist)' * state_dist.final;
    results.var_loggrossy = state_dist.final' * (log(ymat_onxgrid) - results.mean_loggrossy).^2 * yTdist;
    results.var_lognety = state_dist.final' * (log(netymat_onxgrid)- results.mean_lognety).^2 * yTdist;
    
    state_dist.wide = reshape(state_dist.final,[p.nxlong p.nyP p.nyF p.nb]);
    con.final_wide = reshape(con.final,[p.nxlong p.nyP p.nyF p.nb]);
    sav.final_wide = reshape(sav.final,[p.nxlong p.nyP p.nyF p.nb]);
    xgrid.final_wide = reshape(xgrid.final,[p.nxlong p.nyP p.nyF p.nb]);
    sav.orig_wide = reshape(sav.orig,[p.nx p.nyP p.nyF p.nb]);
    con.orig_wide = reshape(con.orig,[p.nx p.nyP p.nyF p.nb]);
    mean_x_check = p.R*(1-p.dieprob)*results.mean_s + results.mean_nety;
    temp = permute(state_dist.wide,[2 1 3 4]);
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
    if abs((meap.nxlongety-results.mean_nety)/meap.nxlongety) > 1e-3
        results.issues = [results.issues,'BadNetIncomeMean'];
    end
    if norm(yPdist_check-yPdist) > 1e-3
        results.issues = [results.issues,'BadGrossIncDist'];
    end
    if min(state_dist.final) < - 0.01
        results.issues = [results.issues,'NegativeStateProbability'];
    end

    %% WEALTH DISTRIBUTION
    temp = sortrows([sav.final state_dist.final]);
    sav.final_sort = temp(:,1);
    state_dist.sort = temp(:,2);
    state_dist.cum = cumsum(state_dist.sort);
    
    results.frac_constrained = (sav.final<=p.borrow_lim)' * state_dist.final;
    results.frac_less5perc_labincome = (sav.final<0.05)' * state_dist.final;
    % wealth percentiles;
    percentiles = [0.1 0.25 0.5 0.9 0.99];
    [state_dist.unique,iu] = unique(state_dist.cum);
    wealthps = interp1(state_dist.unique,sav.final_sort(iu),percentiles,'linear');
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
            idx = (sav.final_sort>=bin) & (sav.final_sort<bin+binwidth);
        else
            idx = (sav.final_sort>=bin) & (sav.final_sort<=bin+binwidth);
        end
        values(ibin) = sum(state_dist.sort(idx));
        ibin = ibin + 1;
    end
    
    
    %% Simulate
    if p.Simulate == 1
        [simulations,ssim] = simulate(p,income,labtaxthresh,sav,con,...
            xgrid,lumptransfer,betacumdist,betacumtrans);
    else
        simulations =[];
    end
    results.simulations = simulations;

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
        plot(xgrid.orig_wide(:,1,iyF,iyb),con.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),con.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',1);
        grid;
        xlim([p.borrow_lim p.xmax]);
        title('Consumption Policy Function');
        legend('Lowest income state','Highest income state');

        % savings policy function
        subplot(2,4,2);
        plot(xgrid.orig_wide(:,1,iyF,iyb),sav.orig_wide(:,1,iyF,iyb)./xgrid.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),sav.orig_wide(:,p.nyP,iyF,iyb)./xgrid.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',1);
        hold on;
        plot(sgrid.short,ones(p.nx,1),'k','LineWidth',0.5);
        hold off;
        grid;
        xlim([p.borrow_lim p.xmax]);
        title('Savings Policy Function s/x');

        % consumption policy function: zoomed in
        subplot(2,4,3);
        plot(xgrid.orig_wide(:,1,iyF),con.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),con.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',2);
        grid;
        xlim([0 4]);
        title('Consumption: Zoomed');

         % savings policy function: zoomed in
        subplot(2,4,4);
        plot(xgrid.orig_wide(:,1,iyF,iyb),sav.orig_wide(:,1,iyF,iyb)./xgrid.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),sav.orig_wide(:,p.nyP,iyF,iyb)./xgrid.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',2);
        hold on;
        plot(sgrid.short,ones(p.nx,1),'k','LineWidth',0.5);
        hold off;
        grid;
        xlim([0 4]);
        title('Savings (s/x): Zoomed');
        
         % gross income distribution
        subplot(2,4,5);
        b = bar(ysort,ysortdist);
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
        title('Asset PMF, Bip.nxlonged');

         % simulation convergence
        if p.Simulate == 1
            subplot(2,4,7);
            plot(1:p.Tsim,mean(ssim),'b','LineWidth',2);
            grid;
            xlim([0 p.Tsim]);
            title('Mean savings (sim)');
        end
    end

    %% COMPUTE SIMULATION MPCs
    if p.ComputeSimMPC ==1
        %theoretical mpc lower bound
        mpclim = p.R*((beta*p.R)^-(1./p.risk_aver))-1;
        Nmpcamount = numel(p.mpcfrac);
        %mpc amounts
        for im = 1:Nmpcamount
            mpcamount{im} = p.mpcfrac{im} * meany;
            xgrid.mpc{im} = xgrid.final_wide + mpcamount{im};
        end

        results.mpcamount = mpcamount;
        
        % Create interpolants for computing mpc's
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
                coninterp{iyP,iyF,ib} = griddedInterpolant(xgrid.final_wide(:,iyP,iyF,ib),con.final_wide(:,iyP,iyF,ib),'linear','linear');
        end
        end
        end
        
        % mpc functions
        for im = 1:Nmpcamount
            mpc{im} = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
            % iterate over (yP,yF,beta)
            for ib = 1:p.nb
            for iyF = 1:p.nyF
            for iyP = 1:p.nyP 
                mpc{im}(:,iyP,iyF,ib) = (coninterp{iyP,iyF,ib}(xgrid.mpc{im}(:,iyP,iyF,ib))...
                    - con.final_wide(:,iyP,iyF,ib))/mpcamount{im};     
            end
            end
            end
            % average mpc, one-period ahead
            results.avg_mpc1{im} = state_dist.final' * mpc{im}(:);
        end
        
        
        
    end
    
    %% COMPUTE MPC FROM STATIONARY DISTRIBUTION
    if p.ComputeDistMPC == 1
        [xgridm,savm,results.avg_mpc,results.DISTmpc_amount] = mpc_forward(xgrid,p,income,sav,...
            betatrans);

    end
    
    %% Print Results
    if p.PrintStats == 1
        print_statistics(results,simulations,p);
    end
end