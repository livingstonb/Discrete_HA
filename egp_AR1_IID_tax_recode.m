% Endogenous Grid Points with AR1 + IID Income
% Cash on Hand as State variable
% Includes NIT and discount factor heterogeneity
% Greg Kaplan 2017

clear;
close all;

path = '/Users/Brian/Documents/GitHub/MPCrecode';
addpath([path '/Auxiliary Functions']);
cd(path);

% Load structure of parameters
prms = load_parameters();
num_prms = size(prms,2);

% Initialize results structure
results = struct();

for ip = 1:num_prms;
    p = prms(ip)
    

    %% INCOME GRIDS

    %persistent income: rowenhurst
    if p.LoadIncomeProcess == 1
        logyPgrid = load('QuarterlyIncomeDynamics/TransitoryContinuous/logyPgrid.txt');
        yPdist = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPdist.txt');
        yPtrans = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPtrans.txt');
        p.nyP = length(logyPgrid);
        logyPgrid = logyPgrid';
    elseif nyP>1
        [logyPgrid, yPtrans, yPdist] = rouwenhorst(p.nyP, -0.5*sd_logyP^2, sd_logyP, rho_logyP);
    else
        logyPgrid = 0;
        yPdist = 1;
        yPtrans = 1;
    end    

    yPgrid = exp(logyPgrid);
    yPcumdist = cumsum(yPdist,1);
    yPcumtrans = cumsum(yPtrans,2);

    % transitory income: disretize normal distribution
    if p.LoadIncomeProcess == 1
        sig2T = load('QuarterlyIncomeDynamics/TransitoryContinuous/sig2T.txt');
        lambdaT = load('QuarterlyIncomeDynamics/TransitoryContinuous/lambdaT.txt');
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

    %     width = fzero(@(x)discrete_normal(nyT,-0.5*sd_logyT^2 ,sd_logyT ,x),2);    
    %     [~,logyTgrid,yTdist] = discrete_normal(nyT,-0.5*sd_logyT^2 ,sd_logyT ,width);

    elseif nyT==1
        logyTgrid = 0;
        yTdist = 1;
    end
    yTgrid = exp(logyTgrid);


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

    % check if yFgrid is row vector like yPgrid, if not, need to rewrite code
    assert(size(yFgrid,1)==1);

    % transition probabilities for yP-yF combined grid
    ytrans = kron(eye(p.nyF),yPtrans);

    p.N = p.nx*p.nyF*p.nyP*p.nb;
    %% DISCOUNT FACTOR

    if p.IterateBeta == 0
        p.maxiterAY = 1;
        beta = p.beta0;
    end

    if p.IterateBeta == 1
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

    % construct matrix of y combinations
    ymat = reshape(repmat(yPgrid',p.ns,1),p.ns*p.nyP,1);
    ymat = repmat(ymat,p.nyF,1) .* reshape(repmat(yFgrid',p.nyP*p.ns,1),p.nyP*p.ns*p.nyF,1);
    ymat = repmat(ymat,p.nb,1)*yTgrid';

    ymatdist = reshape(repmat(yPdist',p.ns,1),p.ns*p.nyP,1);
    ymatdist = repmat(ymatdist,p.nyF,1) .* reshape(repmat(yFdist',p.nyP*p.ns,1),p.nyP*p.ns*p.nyF,1);
    ymatdist = repmat(ymatdist,p.nb,1)*yTdist';

    % find mean y
    ymat_yvals = ymat(1:p.ns:p.ns*p.nyF*p.nyP,:);
    ymatdist_pvals = ymatdist(1:p.ns:p.ns*p.nyF*p.nyP,:);
    temp = sortrows([ymat_yvals(:) ymatdist_pvals(:)],1);
    ysortvals = temp(:,1);
    ycumdist = cumsum(temp(:,2));
    meany = ymat_yvals(:)'*ymatdist_pvals(:);
    totgrossy = meany;

    % find tax threshold on labor income
    if numel(ysortvals)>1
        labtaxthresh = interp1(ycumdist,ysortvals,p.labtaxthreshpc);
    else
        labtaxthresh = 0;
    end    

    totgrossyhigh = max(ymat_yvals(:)-labtaxthresh,0)'*ymatdist_pvals(:);
    lumptransfer = p.labtaxlow*totgrossy + p.labtaxhigh*totgrossyhigh;
    netymat = lumptransfer + (1-p.labtaxlow)*ymat - p.labtaxhigh*max(ymat-labtaxthresh,0);
    netymat_yvals = netymat(1:p.ns:p.ns*p.nyF*p.nyP,:);
    meannety = netymat_yvals(:)'*ymatdist_pvals(:);

    % xgrid, indexed by beta,yF,yP,x
    % cash on hand grid: different min points for each value of (iyP)
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


    % next period's, cash-on-hand as function of saving

    xgrid_wide = reshape(xgrid,p.ns,p.nyP*p.nyF*p.nb);

    if p.IterateBeta == 1
        iterate_EGP = @(x) solve_EGP(x,p,...
    xgrid_wide,ytrans,betatrans,sgrid_wide,u1,u1inv,netymat,meany,...
    yTdist,beq1);

        beta_lb = 1e-5;
        if nb == 1
            beta_ub = betaH - 1e-2;
        else
            beta_ub = betaH - 1e-2 - betawidth;
        end
        results(ip).beta = fmincon(iterate_EGP,beta0,[],[],[],[],beta_lb,beta_ub);
    end



    [AYdiffsq,con,sav,state_dist,cdiff] = solve_EGP(beta,p,...
    xgrid_wide,ytrans,betatrans,sgrid_wide,u1,u1inv,netymat,meany,...
    yTdist,beq1);

    results(ip).mean_s = sav' * state_dist(:);
    results(ip).mean_x = xgrid' * state_dist(:);
    results(ip).mean_y = (ymat*yTdist)' * state_dist(:);
    results(ip).mean_nety = (netymat*yTdist)' * state_dist(:);
    results(ip).mean_x_check = (1+p.r)*results(ip).mean_s + results(ip).mean_nety;
    results(ip).mean_AY = results(ip).mean_s/results(ip).mean_y;
    
    % Record problems
    results(ip).issues = {};
    if cdiff > p.tol_iter
        results(ip).issues = [results(ip).issues,'NoEGPConv'];
    end
    if (results(ip).mean_AY<11) || (results(ip).mean_AY>13)
        results(ip).issues = [results(ip).issues,['BadAY: ' num2str(results(ip).mean_AY)]];
    end
    if abs((results(ip).mean_x-results(ip).mean_x_check)/results(ip).mean_x)> 1e-3
        results(ip).issues = [results(ip).issues,'DistNotStationary'];
    end

    
    results(ip)
    
    %% WEALTH DISTRIBUTION
    results(ip).frac_constrained = (sav<=p.borrow_lim)' * state_dist;

    %% MAKE PLOTS
    newdim = [p.nx p.nyP p.nyF p.nb];
    con_wide = reshape(con,p.nx,p.N/p.nx);
    con_multidim = reshape(con,newdim);
    sav_multidim = reshape(sav,newdim);
    xgrid_multidim = reshape(xgrid,newdim);

    if p.MakePlots ==1 

     figure(1);

        %plot for median fixed effect
        if mod(p.nyF,2)==1
            iyF = (p.nyF+1)/2;
        else
            iyF = iyF/2;
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
        plot(xgrid_multidim(:,1,iyF),con_multidim(:,1,iyF),'b-o',xgrid_multidim(:,p.nyP,iyF),con_multidim(:,p.nyP,iyF),'r-o','LineWidth',2);
        grid;
        xlim([0 4]);
        title('Consumption: Zoomed');

         % savings policy function: zoomed in
        subplot(2,4,4);
        plot(xgrid_multidim(:,1,iyF),sav_multidim(:,1,iyF)./xgrid_multidim(:,1,iyF),'b-o',xgrid_multidim(:,p.nyP,iyF),sav_multidim(:,p.nyP,iyF)./xgrid_multidim(:,p.nyP,iyF),'r-o','LineWidth',2);
        hold on;
        plot(sgrid,ones(p.nx,1),'k','LineWidth',0.5);
        hold off;
        grid;
        xlim([0 4]);
        title('Savings (s/x): Zoomed');
    end

    %% COMPUTE MPCs
    if p.ComputeMPC ==1
        %theoretical mpc lower bound
        mpclim = p.R*((beta*p.R)^-(1./p.risk_aver))-1;
        Nmpcamount = numel(p.mpcfrac);
        %mpc amounts
        for im = 1:Nmpcamount
            mpcamount{im} = p.mpcfrac{im} * meany;
            xgrid_mpc{im} = xgrid + mpcamount{im};
        end
    
        % mpc functions, iterate over (yP,yF,beta)
        for col = 1:p.N/p.nx
            coninterp{col} = griddedInterpolant(xgrid_wide(:,col),con_wide(:,col),'linear');
            mpc{im} = zeros(p.nx,p.N/p.nx);
            for im = 1:Nmpcamount
                mpc{im}(:,col) = (coninterp{col}(xgrid_wide(:,col)+mpcamount{im}) - con_wide(:,col))/mpcamount{im};
            end
        end

        % average mpc
        results(ip).mpcamount = mpcamount;
        for im = 1:Nmpcamount
            results(ip).avg_mpc{im} = mpc{im}(:)' * state_dist(:);
        end
    end
    
    results(ip)
end
