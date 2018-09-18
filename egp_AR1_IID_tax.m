% Endogenous Grid Points with AR1 + IID Income
% Cash on Hand as State variable
% Includes NIT and discount factor heterogeneity
% Greg Kaplan 2017

clear;
close all;

cd('/Users/brianlivingston/Documents/GitHub/MPCrecode')

%% PARAMETERS

Optimized = 0;

%returns
r  = 0.005;
R = 1+ r;

%demographics
dieprob     = 1/50;

% preferences
risk_aver   = 1;
beta        = 0.97365; %if not iterating on discount rate
temptation = 0.0;

betaL = 0.90; %guesses if iterating on discount rate;
betaH = 1/(R*(1-dieprob));

%warm glow bequests: bequest_weight = 0 is accidental
bequest_weight = 0; %0.07;
bequest_luxury = 2; %0.01;

% income risk: AR(1) + IID in logs
LoadIncomeProcess = 0;
nyT         = 5; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

%only relevant if LoadIncomeProcess==0
sd_logyT   = sqrt(0.2431);  %0.20; %relevant if nyT>1
lambdaT = 1; % arrival rate of shocks;

nyP         = 9; %11 persistent component
sd_logyP    = sqrt(0.0047); %0.1950;
rho_logyP   = 0.9947;

nyF     = 1;
sd_logyF    = 0;

% cash on hand / savings grid
nx          = 50;
xmax        = 40;  %multiple of mean gross labor income
xgrid_par   = 0.4; %1 for linear, 0 for L-shaped
borrow_lim  = 0;

%government
labtaxlow       = 0.0; %proportional tax
labtaxhigh      = 0; %additional tax on incomes above threshold
labtaxthreshpc   = 0.99; %percentile of earnings distribution where high tax rate kicks in

savtax          = 0; %0.0001;  %tax rate on savings
savtaxthresh    = 0; %multiple of mean gross labor income


%discount factor shocks;
nb = 2;  %1 or 2
betawidth = 0.065; % beta +/- beta width
betaswitch = 1/50; %0;

% computation
max_iter    = 1000;
tol_iter    = 1.0e-6;
Nsim        = 100000;
Tsim        = 200;

targetAY = 12.0;
maxiterAY = 20;
tolAY = 1.0e-4;

%mpc options
mpcfrac{1}  = 1.0e-10; %approximate thoeretical mpc
mpcfrac{2}  = 0.01; % 1 percent of average gross labor income: approx $500
mpcfrac{3}  = 0.10; % 5 percent of average gross labor income: approx $5000

%% OPTIONS
IterateBeta = 0;
Display     = 1;
MakePlots   = 1;
ComputeMPC  = 1;
SolveDeterministic = 1;

% which function to interpolation 
InterpMUC = 0;
InterpCon = 1; %sometimes useful to stop extrapolating to negative MUC


%% DRAW RANDOM NUMBERS
rng(2017);
yPrand = rand(Nsim,Tsim);
yTrand = randn(Nsim,Tsim);
yFrand = randn(Nsim,1);
betarand = rand(Nsim,Tsim);
dierand = rand(Nsim,Tsim);

%% INCOME GRIDS

%persistent income: rowenhurst
if LoadIncomeProcess == 1
    logyPgrid = load('QuarterlyIncomeDynamics/TransitoryContinuous/logyPgrid.txt');
    yPdist = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPdist.txt');
    yPtrans = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPtrans.txt');
    nyP = length(logyPgrid);
    
elseif nyP>1
    [logyPgrid, yPtrans, yPdist] = rouwenhorst(nyP, -0.5*sd_logyP^2, sd_logyP, rho_logyP);
else
    logyPgrid = 0;
    yPdist = 1;
    yPtrans = 1;
end    
logyPgrid = logyPgrid';
yPgrid = exp(logyPgrid);
yPcumdist = cumsum(yPdist,1);
yPcumtrans = cumsum(yPtrans,2);
    

% transitory income: disretize normal distribution
% if nyT>1
%     width = fzero(@(x)discrete_normal(nyT,-0.5*sd_logyT^2 ,sd_logyT ,x),2);
%     [~,logyTgrid,yTdist] = discrete_normal(nyT,-0.5*sd_logyT^2 ,sd_logyT ,width);
% elseif nyT==1
%     logyTgrid = 0;
%     yTdist = 1;
% end
% yTgrid = exp(logyTgrid);


% transitory income: disretize normal distribution
if LoadIncomeProcess == 1
    sig2T = load('QuarterlyIncomeDynamics/TransitoryContinuous/sig2T.txt');
    lambdaT = load('QuarterlyIncomeDynamics/TransitoryContinuous/lambdaT.txt');
end

if nyT>1
    
    %moments of mixture distribution
    lmu2 = lambdaT.*sd_logyT^2;
    lmu4 = 3.*lambdaT.*(sd_logyT^4);

    %fit thjose moments
    optionsNLLS = optimoptions(@lsqnonlin,'Display','Off');
    lpar = lsqnonlin(@(lp)discretize_normal_var_kurt(lp,nyT,lmu2,lmu4),[2 0.1],[],[],optionsNLLS);
    [lf,lx,lp] = discretize_normal_var_kurt(lpar,nyT,lmu2,lmu4);
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
if nyF>1
    width = fzero(@(x)discrete_normal(nyF,-0.5*sd_logyF^2 ,sd_logyF ,x),2);
    [~,logyFgrid,yFdist] = discrete_normal(nyF,-0.5*sd_logyF^2 ,sd_logyF ,width);
elseif nyF==1
    logyFgrid = 0;
    yFdist = 1;
end
yFgrid = exp(logyFgrid);
yFcumdist = cumsum(yFdist,1);

% Create income grid
if Optimized==1
    % create ymat of dimension N x nyT
    ymat = reshape(repmat(yPgrid,nx,1),nx*nyP,1);
    ymat = repmat(ymat,nyF,1) .* reshape(repmat(yFgrid,nyP*nx,1),nyP*nx*nyF,1);
    ymat = repmat(ymat,nb,1)*yTgrid';
    
    ymatdist = reshape(repmat(yPdist',nx,1),nx*nyP,1);
    ymatdist = repmat(ymatdist,nyF,1) .* reshape(repmat(yFdist',nyP*nx,1),nyP*nx*nyF,1);
    ymatdist = repmat(ymatdist,nb,1)*yTdist';
    
    % grid of unique incomes values
    ymat_yvals = ymat(1:nx:nx*nyF*nyP,:);
    ymatdist_pvals = ymatdist(1:nx:nx*nyF*nyP,:);
    grosslabincgrid = permute(reshape(ymat_yvals,nyP,nyF,nyT),[1 3 2]);
    grosslabincdist = permute(reshape(ymatdist_pvals,nyP,nyF,nyT),[1 3 2]);
else
    %gross labor income grid;
    grosslabincgrid = zeros(nyP,nyT,nyF);
    grosslabincdist= zeros(nyP,nyT,nyF);
    for iyF= 1:nyF
        for iyT= 1:nyT
            for iyP = 1:nyP
                grosslabincgrid(iyP,iyT,iyF) = yPgrid(iyP)*yTgrid(iyT)*yFgrid(iyF);
                grosslabincdist(iyP,iyT,iyF) = yPdist(iyP)*yTdist(iyT)*yFdist(iyF);
            end
        end
    end 
    grosslabincgridvec = grosslabincgrid(:);
    grosslabincdistvec = grosslabincdist(:);
    temp = sortrows([grosslabincgrid(:) grosslabincdist(:)],1);
end

%gross labor income distribution;
temp = sortrows([grosslabincgrid(:) grosslabincdist(:)],1);
glincgridvec = temp(:,1);
glincdistvec = temp(:,2);
glincdistcumvec = cumsum(glincdistvec);


% labor taxes and transfers function
meangrosslabinc = sum(glincgridvec.*glincdistvec);
totgrosslabinc = meangrosslabinc;
if numel(glincgridvec)>1
    labtaxthresh = lininterp1(glincdistcumvec,glincgridvec,labtaxthreshpc);
else
    labtaxthresh = 0;
end    
totgrosslabinchigh = sum(max(glincgridvec - labtaxthresh,0).*glincdistvec);

lumptransfer  = labtaxlow*totgrosslabinc + labtaxhigh*totgrosslabinchigh ;

if Optimized==1
    % net labor income, with dimension N by nyT
    netymat = lumptransfer + (1-labtaxlow)*ymat - labtaxhigh*max(ymat-labtaxthresh,0);
end
    %net labor income
netlabincgrid  = lumptransfer + (1-labtaxlow).*grosslabincgrid - labtaxhigh.*max(grosslabincgrid- labtaxthresh,0);
netlabincgridvec = netlabincgrid(:);
temp = sortrows([netlabincgrid(:) grosslabincdist(:)],1);
nlincgridvec = temp(:,1);
nlincdistvec = temp(:,2);
meannetlabinc = sum(nlincgridvec.*nlincdistvec);


%% SIMULATE INCOME BEFORE SOLVING SINCE EXOGENOUS
if 1
    %initialize 
    yPindsim = zeros(Nsim,Tsim);
    logyTsim = zeros(Nsim,Tsim);
    logyPsim = zeros(Nsim,Tsim);
    logyFsim = zeros(Nsim,1);
    ygrosssim = zeros(Nsim,Tsim);    
    ynetsim = zeros(Nsim,Tsim);    

    %indicators for dying
    diesim = dierand<dieprob;

    %fixed effect;
    [~,idx] = max(yFrand<=yFcumdist',[],2);
    yFindsim = idx;

    %initialize permanent income;
    it = 1;
    [~,idx] = max(yPrand(:,it)<=yPcumdist',[],2);
    yPindsim(:,it) = idx;

    %loop over time periods
    for it = 1:Tsim
        if Display >=1 && mod(it,50) ==0
            disp([' Simulating income realizations, time period ' int2str(it)]);
        end

        %people who dont die: use transition matrix
        if it > 1
                [~,idx] = max(yPrand(diesim(:,it)==0,it)<=yPcumtrans( yPindsim(diesim(:,it)==0,it-1),:),[],2);
                yPindsim(diesim(:,it)==0,it) = idx;
        end

        %people who die: descendendents draw income from stationary distribution
        [~,idx] = max(yPrand(diesim(:,it)==1,it)<=yPcumdist',[],2);
        yPindsim(diesim(:,it)==1,it) = idx;

        logyPsim(:,it) = logyPgrid(yPindsim(:,it));

        %transitory income realization
        if nyT>1
            logyTsim(:,it) = - 0.5*sd_logyT^2 + yTrand(:,it).*sd_logyT;
        else
            logyTsim(:,it) = 0;
        end
        ygrosssim(:,it) = exp(logyPsim(:,it) + logyTsim(:,it)+ logyFsim);
        ynetsim(:,it) = lumptransfer + (1-labtaxlow).*ygrosssim(:,it) - labtaxhigh.*max(ygrosssim(:,it)-labtaxthresh,0);

    end
end
 

%% ASSET GRIDS

% cash on hand grids: different min points for each value of (iyP)
if Optimized==1
    xgrid = linspace(0,1,nx)';
    xgrid = xgrid.^(1/xgrid_par);
    xgrid = repmat(xgrid,nyP*nyF,1);
    xgrid = borrow_lim + (xmax-borrow_lim)*xgrid;
    netymat_minyt = min(netymat(1:nx*nyP*nyF,:),[],2);
    xgrid = xgrid + netymat_minyt;
    xgrid = reshape(xgrid,nx,nyP,nyF);
else
    xgrid = zeros(nx,nyP,nyF);
    for iyF = 1:nyF
    for iyP = 1:nyP
        xgrid(:,iyP,iyF) = linspace(0,1,nx)';
        xgrid(:,iyP,iyF) = xgrid(:,iyP,iyF).^(1./xgrid_par);
        xgrid(:,iyP,iyF) = borrow_lim + min(netlabincgrid(iyP,:,iyF)) + (xmax-borrow_lim).*xgrid(:,iyP,iyF);
    end
    end
end

% savings grid for EGP
ns = nx;
sgrid = linspace(0,1,nx)';
sgrid = sgrid.^(1./xgrid_par);
sgrid = borrow_lim + (xmax-borrow_lim).*sgrid;



%% UTILITY FUNCTION, BEQUEST FUNCTION

if risk_aver==1
    u = @(c)log(c);
    beq = @(a) bequest_weight.* log(a+ bequest_luxury);
else    
    u = @(c)(c.^(1-risk_aver)-1)./(1-risk_aver);
    beq = @(a) bequest_weight.*((a+bequest_luxury).^(1-risk_aver)-1)./(1-risk_aver);
end    

u1 = @(c) c.^(-risk_aver);
u1inv = @(u) u.^(-1./risk_aver);

beq1 = @(a) bequest_weight.*(a+bequest_luxury).^(-risk_aver);


%% INITIALIZE CONSUMPTION FUNCTION

if Optimized==1
    con = r*xgrid;
else
    con = zeros(nx,nyP,nb,nyF);
    for iyF = 1:nyF
    for ib = 1:nb
    for iyP = 1:nyP
        con(:,iyP,ib,iyF) = r.*xgrid(:,iyP,iyF);
    end
    end
    end
end


%% ITERATE ON DISCOUNT FACTOR
if IterateBeta == 0
    maxiterAY = 1;    
end

if IterateBeta == 1
    beta = (betaH + betaL)/2;
end

%initial discount factor grid
if  nb == 1
    betagrid = beta;
    betadist = 1;
    betatrans = 1;
elseif nb ==2 
    betagrid = [beta-betawidth;beta+betawidth];
    betadist = [0.5;0.5];
    betatrans = [1-betaswitch betaswitch; betaswitch 1-betaswitch]; %transitions on average once every 40 years;
else
    error('nb must be 1 or 2');
end
betacumdist = cumsum(betadist);
betacumtrans = cumsum(betatrans,2);

%initialize arrays
mucinterp = cell(nyP,nb,nyF);
coninterp = cell(nyP,nb,nyF);
sav = zeros(nx,nyP,nb,nyF);
muc1next= zeros(ns,nyP,nb,nyT,nyF);
con1next= zeros(ns,nyP,nb,nyT,nyF);
muc1= zeros(ns,nyP,nb,nyF);
con1= zeros(ns,nyP,nb,nyF);
cash1= zeros(ns,nyP,nb,nyF);
temptsc =zeros(nx,nyP,nb,nyF);   

iterAY = 1;
AYdiff = 1;
while iterAY<=maxiterAY && abs(AYdiff)>tolAY
    
    
    
    %% ITERATE ON EULER EQUATION WITH ENDOGENOUS GRID POINTS

    iter = 0;
    cdiff = 1000;

    while iter <= max_iter && cdiff>tol_iter
        iter = iter + 1;

        conlast = con;
        muc  = u1(conlast);   %muc on grid for x' 

        %loop over fixed effect
        for iyF = 1:nyF
            
            %create interpolating function for each yP'
            for ib2 = 1:nb
            for iyP2 = 1:nyP

                if InterpMUC==1
                    mucinterp{iyP2,ib2,iyF} = griddedInterpolant(xgrid(:,iyP2,iyF),muc(:,iyP2,ib2,iyF)- (temptation/(1+temptation)) .* u1(xgrid(:,iyP2,iyF)),'linear');
                    for iyT2 = 1:nyT                
                        muc1next(:,iyP2,ib2,iyT2,iyF) = mucinterp{iyP2,ib2,iyF}(R.*sgrid + netlabincgrid(iyP2,iyT2,iyF));
                    end
                elseif InterpCon==1
                    coninterp{iyP2,ib2,iyF} = griddedInterpolant(xgrid(:,iyP2,iyF),conlast(:,iyP2,ib2,iyF),'linear');
                    for iyT2 = 1:nyT                
                        con1next(:,iyP2,ib2,iyT2,iyF) = coninterp{iyP2,ib2,iyF}(R.*sgrid + netlabincgrid(iyP2,iyT2,iyF));
                        muc1next(:,iyP2,ib2,iyT2,iyF) = u1(con1next(:,iyP2,ib2,iyT2,iyF)) -  (temptation/(1+temptation)) .* u1(R.*sgrid + netlabincgrid(iyP2,iyT2,iyF));
                    end
                end    
            end
            end


            % loop over current persistent income and discount factor
            emuc= zeros(ns,nyP,nb);

            for ib = 1:nb
            for iyP = 1:nyP

                %loop over future income realizations    and discount factor
                for ib2 = 1:nb
                    for iyP2 = 1:nyP
                        for iyT2 = 1:nyT
                                emuc(:,iyP,ib) = emuc(:,iyP,ib) + muc1next(:,iyP2,ib2,iyT2,iyF) * yPtrans(iyP,iyP2) * yTdist(iyT2) * betatrans(ib,ib2);
                        end
                    end
                end
                muc1(:,iyP,ib,iyF) = (1-dieprob) .*  betagrid(ib).*R.*emuc(:,iyP,ib)./(1+savtax.*(sgrid>=savtaxthresh)) + dieprob.*beq1(sgrid);
                con1(:,iyP,ib,iyF) = u1inv(muc1(:,iyP,ib,iyF));
                cash1(:,iyP,ib,iyF) = con1(:,iyP,ib,iyF) + sgrid + savtax.*sgrid.*(sgrid>=savtaxthresh);


                % loop over current period cash on hand
                for ix = 1:nx 
                    if xgrid(ix,iyP,iyF)<cash1(1,iyP,ib,iyF) %borrowing constraint binds
                        sav(ix,iyP,ib,iyF) = borrow_lim;
                    else %borrowing constraint does not bind;
                        sav(ix,iyP,ib,iyF) = lininterp1(cash1(:,iyP,ib,iyF),sgrid,xgrid(ix,iyP,iyF));
                    end                
                end
                con(:,iyP,ib,iyF) = xgrid(:,iyP,iyF) - sav(:,iyP,ib,iyF);
            end   
            end
        end
            cdiff = max(abs(con(:)-conlast(:)));
            if Display >=1 && mod(iter,50) ==0
                disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
            end
    end    


    %% SIMULATE

    xsim = zeros(Nsim,Tsim);
    ssim = zeros(Nsim,Tsim);
    betaindsim = ones(Nsim,Tsim);
    
    %create interpolating functions
    for iyF = 1:nyF
    for ib = 1:nb
    for iyP = 1:nyP
        savinterp{iyP,ib,iyF} = griddedInterpolant(xgrid(:,iyP,iyF),sav(:,iyP,ib,iyF),'linear');
    end
    end
    end
%     fsavinterp = @(iyP,ib,x,iyF)savinterp{iyP,ib,iyF}(x);

    %initialize discount factor;
    if nb>1
        it = 1;
        [~,idx] = max(betarand(:,it)<=betacumdist',[],2);
        betaindsim(:,it) = idx;
    end

    %loop over time periods
    for it = 1:Tsim
        if Display >=1 && mod(it,50) ==0
            disp([' Simulating, time period ' int2str(it)]);
        end

        %discount factor realization
        if nb>1 && it>1
                [~,idx] = max(betarand(:,it)<=betacumtrans( betaindsim(:,it-1),:),[],2);
                betaindsim(:,it) = idx;
        end

        %update cash on hand
        if it>1
            xsim(:,it) = R.*ssim(:,it-1) + ynetsim(:,it);
        end

        %savings choice
%         ssim(:,it) = arrayfun(fsavinterp,yPindsim(:,it),betaindsim(:,it),xsim(:,it));
        for iyF = 1:nyF
        for ib = 1:nb
        for iyP = 1:nyP
            idx = yPindsim(:,it)==iyP & betaindsim(:,it)==ib & yFindsim==iyF; 
            ssim(idx,it) = savinterp{iyP,ib,iyF}(xsim(idx,it));
        end
        end
        end
        ssim(ssim(:,it)<borrow_lim,it) = borrow_lim;
    end

    %convert to assets
    asim = (xsim - ynetsim)./R;

    %consumption
    csim = xsim - ssim - savtax.*max(ssim-savtaxthresh,0);

 
    %% UPDATE BETA ITERATION

       aysim = asim(:,Tsim) ./ meangrosslabinc;
       meanay = mean(aysim);
       
       AYdiff = meanay/targetAY -1;
       if(IterateBeta==1) 
           if AYdiff<0 % too little assets, increase beta
              disp(['Discount factor iteration ' int2str(iterAY) ', beta is ' num2str(beta) ' AY too low: ' num2str(AYdiff)]);
              betaL = beta;
           elseif AYdiff>0 % too much assets, decrease beta
              disp(['Discount factor iteration ' int2str(iterAY) ', beta is ' num2str(beta) ' AY too high: ' num2str(AYdiff)]);
              betaH = beta;
           end
       end
            
       %update discount factor grid
       if iterAY< maxiterAY && abs(AYdiff)>=tolAY
           beta = 0.5*(betaL + betaH);
            if  nb == 1
                betagrid = beta;
                betadist = 1;
                betatrans = 1;
            elseif nb ==2 
                betagrid = [beta-betawidth;beta+betawidth];
                betadist = [0.5;0.5];
                betatrans = [1-betaswitch betaswitch; betaswitch 1-betaswitch]; %transitions on average once every 40 years;
            else
                error('nb must be 1 or 2');
            end
            betacumdist = cumsum(betadist);
            betacumtrans = cumsum(betatrans,2);
       end
       
        iterAY = iterAY+1;
end

%% SOLVE VERSION WITHOUT INCOME RISK FOR COMPARISON
if SolveDeterministic == 1
    
    %cash on hand grid
    norisk.xgrid = linspace(0,1,nx)';
    norisk.xgrid = norisk.xgrid.^(1./xgrid_par);
    norisk.xgrid = borrow_lim + meannetlabinc + (xmax-borrow_lim).*norisk.xgrid;

    
    iter = 0;
    cdiff = 1000;
    
    %initialize arrays
    norisk.mucinterp = cell(nb);
    norisk.coninterp = cell(nb);
    norisk.sav = zeros(nx,nb);
    norisk.muc1next= zeros(ns,nb);
    norisk.muc1= zeros(ns,nb);
    norisk.con1= zeros(ns,nb);
    norisk.cash1= zeros(ns,nb);
   
    %  guess
    for ib  = 1:nb
        norisk.con(:,ib) = r.*norisk.xgrid;
    end
    
    while iter <= max_iter && cdiff>tol_iter
        iter = iter + 1;

        norisk.conlast = norisk.con;
        norisk.muc  = u1(norisk.conlast);   %muc on grid for x' 

        %create interpolating function for each b'
        for ib2 = 1:nb
            if InterpMUC==1
                norisk.mucinterp{ib2} = griddedInterpolant(norisk.xgrid,norisk.muc(:,ib2),'linear');
                norisk.muc1next(:,ib2) = norisk.mucinterp{ib2}(R.*sgrid + meannetlabinc);
                
            elseif InterpCon==1
                norisk.coninterp{ib2} = griddedInterpolant(norisk.xgrid,norisk.conlast(:,ib2),'linear');
                norisk.muc1next(:,ib2) = u1(norisk.coninterp{ib2}(R.*sgrid + meannetlabinc));
            end    
        end


        % loop over discount factor
        norisk.emuc= zeros(ns,nb);

        for ib = 1:nb

            %loop over future discount factor
            for ib2 = 1:nb
                norisk.emuc(:,ib) = norisk.emuc(:,ib) + norisk.muc1next(:,ib2) * betatrans(ib,ib2);
            end
            norisk.muc1(:,ib) = (1-dieprob) .*  betagrid(ib).*R.*norisk.emuc(:,ib)./(1+savtax.*(sgrid>=savtaxthresh)) + dieprob.*beq1(sgrid);
            norisk.con1(:,ib) = u1inv(norisk.muc1(:,ib));
            norisk.cash1(:,ib) = norisk.con1(:,ib) + sgrid + savtax.*sgrid.*(sgrid>=savtaxthresh);


            % loop over current period cash on hand
            for ix = 1:nx 
                if norisk.xgrid(ix)<norisk.cash1(1,ib) %borrowing constraint binds
                    norisk.sav(ix,ib) = borrow_lim;
                else %borrowing constraint does not bind;
                    norisk.sav(ix,ib) = lininterp1(norisk.cash1(:,ib),sgrid,norisk.xgrid(ix));
                end                
            end
            norisk.con(:,ib) = norisk.xgrid - norisk.sav(:,ib);
        end   
     

        cdiff = max(abs(norisk.con(:)-norisk.conlast(:)));
        if Display >=1 && mod(iter,50) ==0
            disp([' EGP Iteration (no risk) ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end    

 end

%% MAKE PLOTS
if MakePlots ==1 
    figure(1);
    
    %plot for median fixed effect
    if mod(nyF,2)==1
        iyF = (nyF+1)/2;
    else
        iyF = iyF/2;
    end
    
    % consumption policy function
    subplot(2,4,1);
    plot(xgrid(:,1,iyF),con(:,1,iyF),'b-',xgrid(:,nyP,iyF),con(:,nyP,iyF),'r-','LineWidth',1);
    grid;
    xlim([borrow_lim xmax]);
    title('Consumption Policy Function');
    legend('Lowest income state','Highest income state');

    % savings policy function
    subplot(2,4,2);
    plot(xgrid(:,1,iyF),sav(:,1,iyF)./xgrid(:,1,iyF),'b-',xgrid(:,nyP,iyF),sav(:,nyP,iyF)./xgrid(:,nyP,iyF),'r-','LineWidth',1);
    hold on;
    plot(sgrid,ones(nx,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([borrow_lim xmax]);
    title('Savings Policy Function s/x');
    
    % consumption policy function: zoomed in
    subplot(2,4,3);
    plot(xgrid(:,1,iyF),con(:,1,iyF),'b-o',xgrid(:,nyP,iyF),con(:,nyP,iyF),'r-o','LineWidth',2);
    grid;
    xlim([0 4]);
    title('Consumption: Zoomed');
    
     % savings policy function: zoomed in
    subplot(2,4,4);
    plot(xgrid(:,1,iyF),sav(:,1,iyF)./xgrid(:,1,iyF),'b-o',xgrid(:,nyP,iyF),sav(:,nyP,iyF)./xgrid(:,nyP,iyF),'r-o','LineWidth',2);
    hold on;
    plot(sgrid,ones(nx,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([0 4]);
    title('Savings (s/x): Zoomed');
    
    
    %income distribution
    subplot(2,4,5);
    hist(ygrosssim(:,Tsim),50);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0 0.5 0.5],'EdgeColor','blue','LineStyle','-');
    ylabel('')
    title('Gross Income distribution');
    
    %asset distribution
    subplot(2,4,6:7);
    hist(asim(:,Tsim),100);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.7 .7 .7],'EdgeColor','black','LineStyle','-');
    ylabel('')
    title('Asset distribution');

    %convergence check
    subplot(2,4,8);
    plot([1:Tsim]',mean(asim,1),'k-','LineWidth',1.5);
    ylabel('Time Period');
    title('Mean Asset Convergence');
    
   % distribution statistics
    disp(['Variance log gross labinc: ' num2str( var(log(ygrosssim(:,Tsim))))]);
    disp(['Mean assets  (relative to mean income) : ' num2str(mean(aysim))]);
    disp(['Fraction borrowing constrained: ' num2str(sum(aysim==borrow_lim)./Nsim * 100) '%']);
    disp(['Fraction <5% mean lab inc: ' num2str(sum(aysim<=0.05)./Nsim *100) '%']); 
    disp(['10th Percentile: ' num2str(quantile(aysim,.1))]);
    disp(['25th Percentile: ' num2str(quantile(aysim,.25))]);
    disp(['50th Percentile: ' num2str(quantile(aysim,.5))]);
    disp(['90th Percentile: ' num2str(quantile(aysim,.9))]);
    disp(['99th Percentile: ' num2str(quantile(aysim,.99))]);

    %table of parameters;
    tabparam  = zeros(4,1);
    tabparam(1) = beta; %discount factor
    tabparam(2) = risk_aver; %risk aversion
    tabparam(3) = labtaxlow; %labor tax rate
    tabparam(4) = lumptransfer./meangrosslabinc; %relative to mean gross lab inc
    
    %table of statistics
    tabstat = zeros(14,1);
    tabstat(1) = var(log(ygrosssim(:,Tsim))); %variance log gross earnings;
    tabstat(2) = ginicoeff(ygrosssim(:,Tsim)); %gini gross earnings;
    tabstat(3) = var(log(ynetsim(:,Tsim))); %variance log net earnings;
    tabstat(4) = ginicoeff(ynetsim(:,Tsim)); %gini net earnings;
    tabstat(5) = var(log(csim(:,Tsim))); %variance log consumption;
    tabstat(6) = ginicoeff(csim(:,Tsim)); %gini consumption;
    tabstat(7) = mean(aysim); %mean wealth (relative to mean earnings);
    tabstat(8) = quantile(aysim,.5); %median wealth (relative to mean earnings);
    tabstat(9) = ginicoeff(aysim); %wealth gini;
    tabstat(10) = quantile(aysim,.9) ./ quantile(aysim,.5); %90-10 wealth distribution
    tabstat(11) = quantile(aysim,.99) ./ quantile(aysim,.5); %99-50 wealth distribution
    tabstat(12) = sum(aysim<=0)./Nsim; %fraction  less than equal to zero;
    tabstat(13) = sum(aysim<=0.05)./Nsim; %fraction  less than equal to 5% av. annual gross earnings;
    tabstat(14) = sum(aysim(aysim>=quantile(aysim,.9))) ./ sum(aysim); %top 10% wealth share;
    tabstat(15) = sum(aysim(aysim>=quantile(aysim,.99))) ./ sum(aysim); %top 1% wealth share;
    tabstat(16) = sum(aysim(aysim>=quantile(aysim,.999))) ./ sum(aysim); %top 0.1% wealth share;
   
end

%% COMPUTE MPCs
if ComputeMPC ==1
    %theoretical mpc lower bound
    mpclim = R*((beta*R)^-(1./risk_aver))-1;
    Nmpcamount = numel(mpcfrac);
    %mpc amounts
    for im = 1:Nmpcamount
        mpcamount{im} = mpcfrac{im}*meangrosslabinc;
        xsim_mpc{im} = xsim(:,Tsim) +mpcamount{im};

    end
    
    %mpc functions
    for iyF = 1:nyF
        for ib = 1:nb
            for iyP = 1:nyP
                %create interpolating function
                coninterp{iyP,ib,iyF} = griddedInterpolant(xgrid(:,iyP,iyF),con(:,iyP,ib,iyF),'linear');
                for im = 1:Nmpcamount
                    mpc{im}(:,iyP,ib,iyF) = ( coninterp{iyP,ib,iyF}(xgrid(:,iyP,iyF)+mpcamount{im}) - con(:,iyP,ib,iyF) ) ./ mpcamount{im};
                end
            end
        end
        %no risk MPCs
        if SolveDeterministic==1
            norisk.coninterp{ib} = griddedInterpolant(norisk.xgrid,norisk.con(:,ib),'linear');
            for im = 1:Nmpcamount
                norisk.mpc{im}(:,ib) = ( norisk.coninterp{ib}(norisk.xgrid+mpcamount{im}) - norisk.con(:,ib) ) ./ mpcamount{im};
            end
        end
    end
    
    % 1 period mpc simulations;
    for im = 1:Nmpcamount
        mpc1sim{im} = zeros(Nsim,1);
        concheck1{im} = zeros(Nsim,1);
        concheck2{im} = zeros(Nsim,1);
        for iyF = 1:nyF
            for ib = 1:nb
                for iyP = 1:nyP
                    idx = (yPindsim(:,Tsim)==iyP & betaindsim(:,Tsim)==ib & yFindsim==iyF);
                    mpc1sim{im}(idx) = ( coninterp{iyP,ib,iyF}(xsim_mpc{im}(idx)) - coninterp{iyP,ib,iyF}(xsim(idx,Tsim)) ) ./ mpcamount{im};
                    concheck1{im}(idx) = coninterp{iyP,ib,iyF}(xsim_mpc{im}(idx));
                    concheck2{im}(idx) = coninterp{iyP,ib,iyF}(xsim(idx,Tsim));
                end
            end
        end
    end
    
    % 4 period mpc simulations;
    
    cbaseline4 = sum(csim(:,Tsim-3:Tsim),2);

    for im = 1:Nmpcamount
        mpc4sim{im} = zeros(Nsim,1);
        
        xtreatsim{im} = zeros(Nsim,4);
        streatsim{im} = zeros(Nsim,4);
        xtreatsim{im}(:,1) = xsim(:,Tsim-3) +mpcamount{im};
    
        for it = 1:4
            %update cash on hand
            if it>1
                xtreatsim{im}(:,it) = R.*streatsim{im}(:,it-1) + ynetsim(:,Tsim-4+it);
            end

            %savings choice
            for iyF = 1:nyF
            for ib = 1:nb
            for iyP = 1:nyP
                idx = (yPindsim(:,Tsim-4+it)==iyP & betaindsim(:,Tsim-4+it)==ib & yFindsim==iyF);
                streatsim{im}(idx,it) = savinterp{iyP,ib,iyF}(xtreatsim{im}(idx,it));
            end
            end
            end
            streatsim{im}(streatsim{im}(:,it)<borrow_lim,it) = borrow_lim;
        end

        ctreatsim{im} = xtreatsim{im} - streatsim{im} - savtax.*max(streatsim{im}-savtaxthresh,0);
        mpc4sim{im} = (sum(ctreatsim{im},2)-cbaseline4)./mpcamount{im};
    end    
    

    
    %wealth distribution and mean MPC by wealth;
%     ayhistgrid = [linspace(borrow_lim,quantile(aysim,.99),100)'; max(aysim)];
    ayhistgrid = [borrow_lim+0.001; (0.01:0.01:0.05)'; (0.055:0.05:0.5)'; (0.6:0.1:1)'; (1.2:0.2:10)'; (11:1:20)'; max(aysim)+0.01];
    aycumdist = sum(aysim<=ayhistgrid',1)'./Nsim;
    aydist = [aycumdist(1); aycumdist(2:end)-aycumdist(1:end-1)];
    
    ayhistgrid0 = [borrow_lim; ayhistgrid(1:end-1)];
    aydelta = ayhistgrid-ayhistgrid0;
%     aydelta = [0.05; ayhistgrid(2:end)-ayhistgrid(1:end-1)];
%     [~,ayidx]= max(aysim<=ayhistgrid' & aysim>ayhistgrid0',[],2);
    [~,ayidx]= max(aysim<ayhistgrid' & aysim>=ayhistgrid0',[],2);
    for im = 1:Nmpcamount
        mpc1dist{im}= accumarray(ayidx,mpc1sim{im})./accumarray(ayidx,1);
        mpc4dist{im}= accumarray(ayidx,mpc4sim{im})./accumarray(ayidx,1);
    end
    
    %MPC in no risk economy on ayhistgrid
    if SolveDeterministic==1
        norisk.agrid = (norisk.xgrid- meannetlabinc)./R;
        norisk.aygrid = norisk.agrid ./ meangrosslabinc;
        
        for im = 1:Nmpcamount
            %average over beta values;
            norisk.avmpc{im} = norisk.mpc{im}*betadist;
            %interpolate on ayhistgrid;
            norisk.mpcinterp{im} = griddedInterpolant(norisk.aygrid,norisk.avmpc{im},'linear');
            norisk.mpc1dist{im} = norisk.mpcinterp{im}(ayhistgrid);
        end                  
    end  
    
    % make MPC figure
 figure;
    [fbar]   = bar(ayhistgrid(1:end-1), [aydist(1); aydist(2:end-1)./aydelta(2:end-1)],'histc');
    axbar = gca;
    sh  = findall(gcf,'marker','*'); delete(sh);
    set(fbar,'FaceColor','blue','EdgeColor','none');
    alpha(0.15);
    axbar.YAxisLocation = 'right';
    if ayhistgrid(end-1)>borrow_lim %to deal with no risk case
%         xlim([borrow_lim ayhistgrid(end-1)]);
        xlim([borrow_lim 5]);
    else
        xlim([borrow_lim borrow_lim+1]);
    end
    set(gca,'FontSize',16) ;
    
    axline = axes('Position',axbar.Position,'XAxisLocation','bottom','YAxisLocation','left','Color','none');
    hold on;
    fmpc1 = plot(ayhistgrid(1:end-1),mpc1dist{1}(1:end-1),'Parent',axline);
    set(fmpc1,'LineWidth',2,'Color',[0.1 0.1 0.9],'LineStyle',':');
    fmpc2 = plot(ayhistgrid(1:end-1),mpc1dist{2}(1:end-1),'Parent',axline);
    set(fmpc2,'LineWidth',2,'Color',[0.9 0.1 0.1]);    
    fmpc3 = plot(ayhistgrid(1:end-1),mpc1dist{3}(1:end-1),'Parent',axline);
    set(fmpc3,'LineWidth',2,'Color',[0.1 0.9 0.1],'LineStyle','--');

    fmpc4per2 = plot(ayhistgrid(1:end-1),mpc4dist{2}(1:end-1),'Parent',axline);
    set(fmpc4per2,'LineWidth',2,'Color',[0.4 0.4 0.4],'LineStyle','-.');

    fmpcnorisk = plot(ayhistgrid(1:end-1),norisk.mpc1dist{2}(1:end-1),'Parent',axline);
    set(fmpcnorisk,'LineWidth',1,'Color',[0.6 0 0.6],'LineStyle','-');

    lgd = legend({['1 per. MPC, theoretical: mean= ' num2str(round(100*sum(aydist.*mpc1dist{1}),1)) '%'],...
                            ['1 per. MPC, ' num2str(mpcfrac{2}) 'x av. lab inc: mean=  ' num2str(round(100*sum(aydist.*mpc1dist{2}),1)) '%'],...
                            ['1 per. MPC, ' num2str(mpcfrac{3}) 'x av. lab inc: mean= ' num2str(round(100*sum(aydist.*mpc1dist{3}),1)) '%'],...
                            ['4 per. MPC, ' num2str(mpcfrac{2}) 'x av. lab inc: mean= ' num2str(round(100*sum(aydist.*mpc4dist{2}),1)) '%'],...
                            ['1 per. MPC, ' num2str(mpcfrac{2}) 'x av. lab inc (no risk): mean= ' num2str(round(100*sum(aydist.*norisk.mpc1dist{2}),1)) '%'],...
                });
    
    grid on;
    hold off;
    if ayhistgrid(end-1)>borrow_lim %to deal with no risk case
%         xlim([borrow_lim ayhistgrid(end-1)]);
        xlim([borrow_lim 5]);
    else
        xlim([borrow_lim borrow_lim+1]);
    end
set(gca,'FontSize',16) ;
    
    
    
    % mpc distribution statistics
    disp(['Mean 1 period MPC amount 1: ' num2str(mean(mpc1sim{1}))]);
    disp(['Mean 1 period MPC amount 2: ' num2str(mean(mpc1sim{2}))]);
    disp(['Mean 1 period MPC amount 3: ' num2str(mean(mpc1sim{3}))]);
    disp(['Mean 4 period MPC amount 1: ' num2str(mean(mpc4sim{1}))]);
    disp(['Mean 4 period MPC amount 2: ' num2str(mean(mpc4sim{2}))]);
    disp(['Mean 4 period MPC amount 3: ' num2str(mean(mpc4sim{3}))]);
    
    %mpc decomposition
    mpc_true = zeros(1,Nmpcamount);
    mpc_decomp = zeros(4,Nmpcamount);
    for im = 1:Nmpcamount
        mpc_true(1,im) = sum(aydist.*mpc1dist{im});
        mpc_decomp(1,im)= mpclim;
%     mpc_decomp2= sum(aydist.*(norisk.mpc1dist{1}-mpclim));
%     mpc_decomp3= sum(aydist.*(mpc1dist{1}-norisk.mpc1dist{1}));
        mpc_decomp(2,im) = aydist(1).*norisk.mpc1dist{im}(1)-mpclim;
        mpc_decomp(3,im) = sum(aydist(2:end).*norisk.mpc1dist{im}(2:end));
        mpc_decomp(4,im) = sum(aydist.*(mpc1dist{im}-norisk.mpc1dist{im}));
    end
    disp(['MPC theoretical decomp, RA: ' num2str(round(100*mpc_decomp(1,1)./mpc_true(1),1)) '%']);
    disp(['MPC theoretical decomp, borrowing constraints: ' num2str(round(100*mpc_decomp(2,1)./mpc_true(1),1)) '%']);
    disp(['MPC theoretical decomp, income risk (wealth dist): ' num2str(round(100*mpc_decomp(3,1)./mpc_true(1),1)) '%']);
    disp(['MPC theoretical decomp, income risk (mpc function): ' num2str(round(100*mpc_decomp(4,1)./mpc_true(1),1)) '%']);
    
    %compare with SCF dist
    scfdist = load('scf2016_networth.txt');
    for im = 1:Nmpcamount
        mpc_scfdist{im} = sum(scfdist(:,2).*mpc1dist{im});
        mpc_norisk_scfdist{im} = sum(scfdist(:,2).*norisk.mpc1dist{im});
    end
    disp(['Mean 1 period MPC amount 1, SCF dist: ' num2str(mean(mpc_scfdist{1}))]);
    disp(['Mean 1 period MPC amount 2, SCF dist: ' num2str(mean(mpc_scfdist{2}))]);
    disp(['Mean 1 period MPC amount 3, SCF dist: ' num2str(mean(mpc_scfdist{3}))]);
    
    %table of mpc stats;
    tabmpc  = zeros(4,1);
    tabmpc(1) = mean(mpc1sim{1});
    tabmpc(2) = mean(mpc1sim{2});
    tabmpc(3) = mean(mpc1sim{3});
    
end    
    
    
 %% MAKE TABLE
 
  
    
    taboutput = [tabparam; tabstat];
    
    
    
  %% FUNCTIONS

function [f,lx,lp] = discretize_normal_var_kurt(y,n,mu2,mu4)
    lwidth = y(1);
    llambda = y(2);
    lx = linspace(-lwidth*sqrt(mu2), lwidth*sqrt(mu2),n)';
    lp = discrete_normal_alt(lx,0,sqrt(mu2));
    lmass0 = zeros(n,1);
    lmass0((n+1)/2) = 1; 
    
    lp = (1-llambda).*lp + lmass0.*llambda;
    
    Ex2 = sum(lp.*(lx.^2));
    Ex4 = sum(lp.*(lx.^4));
    
    f = [sqrt(Ex2)-sqrt(mu2); Ex4.^0.25 - mu4.^0.25];
    
end    

function [grid, trans, dist] = rouwenhorst(n, mu, sigma, rho)
    
    % grid
    width = sqrt((n-1) * sigma^2 / ( 1 - rho^2));
    grid = linspace( mu-width, mu + width, n)';
    
    %transition matrix
    p0 = (1 + rho) / 2;
    trans = [p0 1-p0; 1-p0 p0];
    
    if n > 2
        for i = 1:n-2
            cstr_temp = zeros(length(trans(:,1)), 1);
            trans = p0 .* [trans cstr_temp; cstr_temp.' 0] + (1 - p0 ) .* [cstr_temp trans; cstr_temp.' 0]  + (1 - p0 ) .*  [ cstr_temp.' 0; trans cstr_temp] + p0 .* [ cstr_temp.' 0; cstr_temp trans];
        end
        trans = bsxfun(@rdivide, trans, sum(trans,2));
    end
    
    % ergodic distribution
    dist = ones(1,n)./n;
    for i = 1: 100
        dist = dist*(trans^i);
    end
    dist = dist';
end
    
    