clear

path = '/Users/Brian/Documents/GitHub/MPCrecode';
addpath([path '/Auxiliary Functions']);

%% SET PARAMETERS
Nsim = 1e6;

% Baseline quarterly
params(1).IncomeProcess = '';
params(1).freq         = 4;
params(1).nyT          = 11;
params(1).yTContinuous = 0;
params(1).sd_logyT     = sqrt(0.2087);
params(1).lambdaT      = 1;
params(1).nyP          = 11;
params(1).sd_logyP     = sqrt(0.0108);
params(1).rho_logyP    = 0.9881;
params(1).nyF          = 1;
params(1).sd_logyF     = 0;
params(1).labtaxlow       = 0; %proportional tax
params(1).labtaxhigh      = 0; %additional tax on incomes above threshold
params(1).labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in
params(1).savtax          = 0; %0.0001;  %tax rate on savings
params(1).savtaxthresh    = 0; %multiple of mean gross labor income
income{1} = gen_income_variables(params(1));

% quarterly_a
params(2) = params(1);
params(2).IncomeProcess = 'IncomeGrids/quarterly_a.mat';
income{2} = gen_income_variables(params(2));

% quarterly_b
params(3) = params(1);
params(3).IncomeProcess = 'IncomeGrids/quarterly_b.mat';
income{3} = gen_income_variables(params(3));

% quarterly_b
params(4) = params(1);
params(4).IncomeProcess = 'IncomeGrids/quarterly_c.mat';
income{4} = gen_income_variables(params(4));

%% SIMULATE
yTsim = cell(1,4);
yPsim = cell(1,4);

for i = 1:4
    yTrand = rand(Nsim,4);
    yPrand = rand(Nsim,4);
    yTindsim = zeros(Nsim,4,'int8');
    yPindsim = zeros(Nsim,4,'int8');
    
    for it = 1:4
        [~,yTindsim(:,it)] = max(bsxfun(@lt,yTrand(:,it),income{i}.yTcumdist'),[],2);
        if it == 1
            [~,yPindsim(:,it)] = max(bsxfun(@lt,yPrand(:,it),income{i}.yPcumdist'),[],2);
        else
            [~,yPindsim(:,it)] = max(bsxfun(@lt,yPrand(:,it),income{i}.yPcumtrans(yPindsim(:,it-1),:)),[],2);
        end
    end
    
    yTsim{i} = income{i}.yTgrid(yTindsim);
    yPsim{i} = income{i}.yPgrid(yPindsim);
end

%% COMPUTE MOMENTS

Qmoments(1).name = 'q_base';
Qmoments(2).name = 'q_a';
Qmoments(3).name = 'q_b';
Qmoments(4).name = 'q_c';
Amoments = Qmoments;

% quarterly moments
for i = 1:4
    yTsimQ = yTsim{i}(:,1);
    Qmoments(i).yT_mu1 = mean(yTsimQ);
    Qmoments(i).yT_mu2 = mean( (yTsimQ - mean(yTsimQ)).^2 );
    Qmoments(i).yT_mu3 = mean( (yTsimQ - mean(yTsimQ)).^3 );
    Qmoments(i).yT_mu4 = mean( (yTsimQ - mean(yTsimQ)).^4 );
    
    logyTsimQ = log(yTsimQ);
    Qmoments(i).logyT_mu1 = mean(logyTsimQ);
    Qmoments(i).logyT_mu2 = mean( (logyTsimQ - mean(logyTsimQ)).^2 );
    Qmoments(i).logyT_mu3 = mean( (logyTsimQ - mean(logyTsimQ)).^3 );
    Qmoments(i).logyT_mu4 = mean( (logyTsimQ - mean(logyTsimQ)).^4 );
    
    yPsimQ = yPsim{i}(:,1);
    Qmoments(i).yP_mu1 = mean(yPsimQ);
    Qmoments(i).yP_mu2 = mean( (yPsimQ - mean(yPsimQ)).^2 );
    Qmoments(i).yP_mu3 = mean( (yPsimQ - mean(yPsimQ)).^3 );
    Qmoments(i).yP_mu4 = mean( (yPsimQ - mean(yPsimQ)).^4 );
    
    logyPsimQ = log(yPsimQ);
    Qmoments(i).logyP_mu1 = mean(logyPsimQ);
    Qmoments(i).logyP_mu2 = mean( (logyPsimQ - mean(logyPsimQ)).^2 );
    Qmoments(i).logyP_mu3 = mean( (logyPsimQ - mean(logyPsimQ)).^3 );
    Qmoments(i).logyP_mu4 = mean( (logyPsimQ - mean(logyPsimQ)).^4 );
end

% annual moments
for i = 1:4
    yTsimA = sum(yTsim{i},2);
    Amoments(i).yT_mu1 = mean(yTsimA);
    Amoments(i).yT_mu2 = mean( (yTsimA - mean(yTsimA)).^2 );
    Amoments(i).yT_mu3 = mean( (yTsimA - mean(yTsimA)).^3 );
    Amoments(i).yT_mu4 = mean( (yTsimA - mean(yTsimA)).^4 );
    
    logyTsimA = log(yTsimA);
    Amoments(i).logyT_mu1 = mean(logyTsimA);
    Amoments(i).logyT_mu2 = mean( (logyTsimA - mean(logyTsimA)).^2 );
    Amoments(i).logyT_mu3 = mean( (logyTsimA - mean(logyTsimA)).^3 );
    Amoments(i).logyT_mu4 = mean( (logyTsimA - mean(logyTsimA)).^4 );
    
    yPsimA = sum(yPsim{i},2);
    Amoments(i).yP_mu1 = mean(yPsimA);
    Amoments(i).yP_mu2 = mean( (yPsimA - mean(yPsimA)).^2 );
    Amoments(i).yP_mu3 = mean( (yPsimA - mean(yPsimA)).^3 );
    Amoments(i).yP_mu4 = mean( (yPsimA - mean(yPsimA)).^4 );
    
    logyPsimA = log(yPsimA);
    Amoments(i).logyP_mu1 = mean(logyPsimA);
    Amoments(i).logyP_mu2 = mean( (logyPsimA - mean(logyPsimA)).^2 );
    Amoments(i).logyP_mu3 = mean( (logyPsimA - mean(logyPsimA)).^3 );
    Amoments(i).logyP_mu4 = mean( (logyPsimA - mean(logyPsimA)).^4 );
end
