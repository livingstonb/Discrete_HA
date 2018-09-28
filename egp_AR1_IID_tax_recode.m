function [simulations,results] = egp_AR1_IID_tax_recode(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % STRUCTURES:
    % basemodel - stores objects from the model with income risk
    % norisk - stores objects from the model without income risk
    % simulations - stores simulation results for basemodel
    % sgrid - stores savings grid in various sizes
    % xgrid - stores cash-on-hand grid in various sizes
    % income - stores objects from the income process
    % prefs - stores objects related to preferences
    
    %% ADJUST PARAMETERS FOR DATA FREQUENCY, OTHER CHOICES

    p.r = p.r/p.freq;
    p.R = 1 + p.r;
    p.beta0 = p.beta0^(1/p.freq);
    p.dieprob = p.dieprob^(1/p.freq);
    p.betaswitch = p.betaswitch^(1/p.freq);
    p.beta = p.betaL^(1/p.freq);
    p.betaH = 1/((p.R)*(1-p.dieprob));
    
    if p.Simulate == 0 && p.ComputeSimMPC == 1
        p.ComputeSimMPC = 0;
        disp('SimMPC turned off because Simulate == 0')
    end
    

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
        results.beta = beta;
        p.beta = beta;
    end

    %initial discount factor grid
    if  p.nb == 1
        prefs.betadist = 1;
        prefs.betatrans = 1;
    elseif p.nb ==2 
        prefs.betadist = [0.5;0.5];
        prefs.betatrans = [1-p.betaswitch p.betaswitch; p.betaswitch 1-p.betaswitch]; %transitions on average once every 40 years;
    else
        error('nb must be 1 or 2');
    end
    prefs.betacumdist = cumsum(prefs.betadist);
    prefs.betacumtrans = cumsum(prefs.betatrans,2);
    

    %% ASSET AND INCOME GRIDS

    sgrid.orig = linspace(0,1,p.nx)';
    sgrid.orig = sgrid.orig.^(1./p.xgrid_par);
    sgrid.orig = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid.orig;
    sgrid.short = sgrid.orig;
    sgrid.norisk_wide = repmat(sgrid.short,[1 p.nb]);
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
    meannety = netymat(:)'*ymatdist(:);

    % xgrid, indexed by beta,yF,yP,x (N by 1 matrix)
    % cash on hand grid: different min points for each value of (iyP,iyF)
    minyT = repmat(kron(min(netymat,[],2),ones(p.nx,1)),p.nb,1);
    xgrid.orig = sgrid.wide(:) + minyT;
    xgrid.orig_wide = reshape(xgrid.orig,[p.nx p.nyP p.nyF p.nb]);
    xgrid.norisk_short = sgrid.short + meannety;
    xgrid.norisk_wide = repmat(xgrid.norisk_short,1,p.nb);
    
    % Store income variables in a structure
    newfields = {'ymat','netymat','meany','original_meany','yPgrid',...
        'yTgrid','yFgrid','yPdist','yTdist','yFdist','yPcumtrans',...
        'yPtrans','yPcumdist','yFcumdist','yTcumdist','ytrans','meannety'};
    for i = 1:numel(newfields)
        income.(newfields{i}) = eval(newfields{i});
    end

    %% UTILITY FUNCTION, BEQUEST FUNCTION

    if p.risk_aver==1
        prefs.u = @(c)log(c);
        prefs.beq = @(a) p.bequest_weight.* log(a+ p.bequest_luxury);
    else    
        prefs.u = @(c)(c.^(1-p.risk_aver)-1)./(1-p.risk_aver);
        prefs.beq = @(a) p.bequest_weight.*((a+p.bequest_luxury).^(1-p.risk_aver)-1)./(1-p.risk_aver);
    end    

    prefs.u1 = @(c) c.^(-p.risk_aver);
    prefs.u1inv = @(u) u.^(-1./p.risk_aver);
    prefs.beq1 = @(a) p.bequest_weight.*(a+p.bequest_luxury).^(-p.risk_aver);


    %% MODEL SOLUTION

    if p.IterateBeta == 1
        [beta,exitflag] = iterate_beta(p,xgrid,sgrid,prefs,income);
        if exitflag ~= 1
            results = struct();
            results.issues = {'NoBetaConv'};
            return
        end
        results.beta = beta;
        p.beta = beta;
    end

    % Use final beta to get policy functions and distribution
    ergodic_tol = 1e-7;
    [~,basemodel,xgrid.longgrid] = solve_EGP(beta,p,xgrid,sgrid,prefs,...
                                            ergodic_tol,income,p.nxlong);
                                        
    xgrid.longgrid_wide = reshape(xgrid.longgrid,[p.nxlong,p.nyP,p.nyF,p.nb]);
    
    %% Store important moments
    
    ymat_onxgrid = repmat(kron(ymat,ones(p.nxlong,1)),p.nb,1);
    netymat_onxgrid = repmat(kron(netymat,ones(p.nxlong,1)),p.nb,1);
    
    results.mean_s = basemodel.sav_longgrid' * basemodel.SSdist;
    results.mean_x = xgrid.longgrid(:)' * basemodel.SSdist;
    results.mean_grossy = (ymat_onxgrid*yTdist)' * basemodel.SSdist;
    results.mean_loggrossy = (log(ymat_onxgrid)*yTdist)' * basemodel.SSdist;
    results.mean_nety = (netymat_onxgrid*yTdist)' * basemodel.SSdist;
    results.mean_lognety = (log(netymat_onxgrid)*yTdist)' * basemodel.SSdist;
    results.var_loggrossy = basemodel.SSdist' * (log(ymat_onxgrid) - results.mean_loggrossy).^2 * yTdist;
    results.var_lognety = basemodel.SSdist' * (log(netymat_onxgrid)- results.mean_lognety).^2 * yTdist;
    
    if p.WealthInherited == 1
        results.mean_x_check = p.R*results.mean_s + results.mean_nety;
    else
        results.mean_x_check = p.R*(1-p.dieprob)*results.mean_s + results.mean_nety;
    end
    temp = permute(basemodel.SSdist_wide,[2 1 3 4]);
    % reconstruct yPdist from computed stationary distribution for error
    % checking
    yPdist_check = sum(sum(sum(temp,4),3),2);
    

    %% Store problems
    results.issues = {};
    if basemodel.EGP_cdiff > p.tol_iter
        results.issues = [results.issues,'NoEGPConv'];
    end
    if  abs((p.targetAY - results.mean_s)/p.targetAY) > 1e-3
        results.issues = [results.issues,'BadAY'];
    end
    if abs((results.mean_x-results.mean_x_check)/results.mean_x)> 1e-3
        results.issues = [results.issues,'DistNotStationary'];
    end
    if abs((income.meannety-results.mean_nety)/income.meannety) > 1e-3
        results.issues = [results.issues,'BadNetIncomeMean'];
    end
    if norm(yPdist_check-yPdist) > 1e-3
        results.issues = [results.issues,'Bad_yP_Dist'];
    end
    if min(basemodel.SSdist) < - 1e-3
        results.issues = [results.issues,'NegativeStateProbability'];
    end

    %% WEALTH DISTRIBUTION
    
    results.frac_constrained = (basemodel.sav_longgrid<=p.borrow_lim)' * basemodel.SSdist;
    results.frac_less5perc_labincome = (basemodel.sav_longgrid<0.05*income.meany)' * basemodel.SSdist;
    % wealth percentiles;
    percentiles = [0.1 0.25 0.5 0.9 0.99];
    [basemodel.SSdist_unique,iu] = unique(basemodel.SScumdist);
    wealthps = interp1(basemodel.SSdist_unique,basemodel.sav_longgrid_sort(iu),percentiles,'linear');
    results.p10wealth = wealthps(1);
    results.p25wealth = wealthps(2);
    results.p50wealth = wealthps(3);
    results.p90wealth = wealthps(4);
    results.p99wealth = wealthps(5);
    
    % create histogram of asset holdings
    binwidth = 0.25;
    bins = 0:binwidth:p.xmax;
    values = zeros(p.xmax+1,1);
    ibin = 1;
    for bin = bins
        if bin < p.xmax
            idx = (basemodel.sav_longgrid_sort>=bin) & (basemodel.sav_longgrid_sort<bin+binwidth);
        else
            idx = (basemodel.sav_longgrid_sort>=bin) & (basemodel.sav_longgrid_sort<=bin+binwidth);
        end
        values(ibin) = sum(basemodel.SSdist_sort(idx));
        ibin = ibin + 1;
    end
    
    %% EGP FOR MODEL WITHOUT INCOME RISK
    % Deterministic model
    if p.SolveDeterministic == 1
        [norisk,norisk.EGP_cdiff] = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income);
    end
    
    %% SIMULATIONS
    % Full model
    if p.Simulate == 1
        [simulations,ssim] = simulate(p,income,labtaxthresh,basemodel,...
                                        xgrid,lumptransfer,prefs,results);
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
        plot(xgrid.orig_wide(:,1,iyF,iyb),basemodel.con_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.con_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',1);
        grid;
        xlim([p.borrow_lim p.xmax]);
        title('Consumption Policy Function');
        legend('Lowest income state','Highest income state');

        % savings policy function
        subplot(2,4,2);
        plot(xgrid.orig_wide(:,1,iyF,iyb),basemodel.sav_wide(:,1,iyF,iyb)./xgrid.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.sav_wide(:,p.nyP,iyF,iyb)./xgrid.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',1);
        hold on;
        plot(sgrid.short,ones(p.nx,1),'k','LineWidth',0.5);
        hold off;
        grid;
        xlim([p.borrow_lim p.xmax]);
        title('Savings Policy Function s/x');

        % consumption policy function: zoomed in
        subplot(2,4,3);
        plot(xgrid.orig_wide(:,1,iyF),basemodel.con_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.con_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',2);
        grid;
        xlim([0 4]);
        title('Consumption: Zoomed');

         % savings policy function: zoomed in
        subplot(2,4,4);
        plot(xgrid.orig_wide(:,1,iyF,iyb),basemodel.sav_wide(:,1,iyF,iyb)./xgrid.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.sav_wide(:,p.nyP,iyF,iyb)./xgrid.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',2);
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
 
    
    %% MPCS
    
    %mpc amounts
    for im = 1:numel(p.mpcfrac)
        mpcamount{im} = p.mpcfrac{im} * income.meany;
        xgrid.mpc{im} = xgrid.longgrid_wide + mpcamount{im};
    end
    
    for im = 1:numel(p.mpcfrac)
        mpc{im} = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        % iterate over (yP,yF,beta)
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP 
            mpc{im}(:,iyP,iyF,ib) = (basemodel.coninterp{iyP,iyF,ib}(xgrid.mpc{im}(:,iyP,iyF,ib))...
                        - basemodel.con_longgrid_wide(:,iyP,iyF,ib))/mpcamount{im};     
        end
        end
        end
        % average mpc, one-period ahead
        results.avg_mpc1{im} = basemodel.SSdist' * mpc{im}(:);
    end
    
        
    for im = 1:numel(p.mpcfrac)
        results.mpcamount{im} = p.mpcfrac{im};
    end
        
    if p.ComputeDistMPC == 1
        [xgridm,savm,results.avg_mpc1_alt,results.avg_mpc4,results.distMPCamount] ...
            = mpc_forward(xgrid,p,income,basemodel,prefs);
    end
    
    %% Print Results
    if p.PrintStats == 1
        print_statistics(results,simulations,p);
    end
    
    results.xgrid = xgrid;
end