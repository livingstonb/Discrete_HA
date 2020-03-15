function [params, all_names] = parameters(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import aux.set_shared_fields

    dollars = [-1, -500, -5000, 1, 500, 5000];
    shared_params.shocks = dollars ./ 72000;

    shared_params.xgrid_par = 0.2;
    shared_params.xgrid_term1wt = 0;
    shared_params.xgrid_term1curv = 1;

    shared_params.shocks_labels = {};
    for ishock = 1:6
        val = dollars(ishock);
        if val < 0
            shared_params.shocks_labels{ishock} = sprintf('-$%g', abs(val));
        else
            shared_params.shocks_labels{ishock} = sprintf('$%g', abs(val));
        end
    end
    % shared_params = struct();

    % location of baseline income process for quarterly case
    quarterly_b_path = 'input/income_quarterly_b_contyT.mat';
    quarterly_c_path = 'input/income_quarterly_c_contyT.mat';

    quarterly_b_params = shared_params;
    quarterly_b_params.sd_logyT = sqrt(0.6376);
    quarterly_b_params.lambdaT = 0.25;

    quarterly_c_params = shared_params;
    quarterly_c_params.sd_logyT = sqrt(1.6243);
    quarterly_c_params.lambdaT = 0.0727;

    quarterly_a_params = shared_params;
    quarterly_a_params.sd_logyT = sqrt(0.2087);
    quarterly_a_params.sd_logyP = sqrt(0.01080);
    quarterly_a_params.rho_logyP = 0.9881;
    quarterly_a_params.lambdaT = 1;

    annual_params = shared_params;
    annual_params.sd_logyT = sqrt(0.0494);
    annual_params.sd_logyP = sqrt(0.0422);
    annual_params.rho_logyP = 0.9525;
    annual_params.lambdaT = 1;
    
    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    
    % Annual
    params(1) = setup.Params(1, 'Annual', '');
    params(1).beta0 = 0.984108034755346;
    params(1) = set_shared_fields(params(1), annual_params);
    params(1).group = {'Baseline', 'A1'};
%     
    % Quarterly
    params(end+1) = setup.Params(4, 'Quarterly', quarterly_b_path);
    params(end) = set_shared_fields(params(end), quarterly_b_params);
    params(end).beta0 = 0.984363510593659;
    params(end).group = {'Baseline', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6'};
    
    %----------------------------------------------------------------------
    % PART 2, DIFFERENT ASSUMPTIONS
    %----------------------------------------------------------------------
    for ifreq = [4]
        if ifreq == 1
            lfreq = 'A';
            IncomeProcess = '';
            income_params = annual_params;
        else
            lfreq = 'Q';
            IncomeProcess = quarterly_b_path;
            income_params = quarterly_b_params;
        end

        % different mean wealth targets
        for mw = [1, 0.5, 0.25]
            name = sprintf('A/Y = %g', mw);
            params(end+1) = setup.Params(4, name, quarterly_b_path);
            params(end) = set_shared_fields(params(end), quarterly_b_params);
            params(end).target_value = mw;
            params(end).group = {'Q1'};
            params(end).betaL = 0.5;
        end

        % no death
        name = 'No Death';
        params(end+1) = setup.Params(4, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end).group = {'Q1'};
        params(end).dieprob = 0;
        params(end).beta0 = 0.975363510593659;

        % no bequests
        name = 'No Bequests';
        params(end+1) = setup.Params(4, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end).group = {'Q1'};
        params(end).Bequests = 0;

        % perfect annuities
        name = 'Annuities';
        params(end+1) = setup.Params(ifreq, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end).group = {'Q1'};
        params(end).annuities = true;
        params(end).betaH0 = - 5e-3;

        % different interest rates
        for ii = [0, 5]
            name = sprintf('r = %g%% p.a.', ii);
            params(end+1) = setup.Params(ifreq, name, quarterly_b_path);
            params(end) = set_shared_fields(params(end), quarterly_b_params);
            params(end).r = ii/100;
            params(end).group = {'Q6'};
            
            if ii == 0
                params(end).label = {'Low r'};
            else
                params(end).label = {'High r'};
                params(end).betaH0 = -3e-3;
                params(end).beta0 = 0.8;
            end
        end
        
        % interest rate heterogeneity
        name = 'Permanent r het, r in {0,2,4} p.a.';
        params(end+1) = setup.Params(ifreq, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end).r = [0, 2, 4] / 100;
        params(end).betaH0 = -1e-4;
        params(end).beta0 = 0.973149481985717;
        params(end).group = {'Q6'};
        params(end).label = {'Heterogeneity in r'};
        params(end).other = {'{0, 2, 4}'};
        
        name = 'Permanent r het, r in {-2,2,6} p.a.';
        params(end+1) = setup.Params(ifreq,name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end).r = [-2, 2, 6] / 100;
        params(end).betaH0 = -1e-4;
        params(end).beta0 = 0.955885729527277;
        params(end).group = {'Q6'};
        params(end).label = {'Heterogeneity in r'};
        params(end).other = {'{-2,2,6}'};


%         % different tax rates
%         for itax = [0.05, 0.1, 0.15, 0.25]
%             name = [lfreq ' LabTax' num2str(itax)];
%             params(end+1) = setup.Params(ifreq,name,'');
%             params(end).labtaxlow = itax;
%         end

        

%         % bequest curvature
%         for bcurv = [0.1 0.5 1 2 5]
%             name = [lfreq ' BeqWt0.02 BeqLux0.01 BeqCurv' num2str(bcurv)];
%             params(end+1) = setup.Params(ifreq,name,IncomeProcess);
%             params(end).bequest_weight = 0.02;
%             params(end).bequest_luxury = 0.01;
%             params(end).bequest_curv   = bcurv;
%         endR
       
        for deathp = 1/50 %[0 1/50]
            if deathp == 0
                deathind = ' NoDeath';
            else
                deathind = ' Death';
            end
             % fixed beta heterogeneity
            for ibw = [0.001, 0.005, 0.01]
                name = sprintf('Beta5, pSwitch0, pSpacing%g', ibw);
                params(end+1) = setup.Params(ifreq, name, IncomeProcess);
                params(end) = set_shared_fields(params(end), income_params);
                params(end).nbeta = 5;
                params(end).betawidth = ibw;
                params(end).prob_zswitch = 0;
                params(end).dieprob = deathp;
                params(end).beta0 = 0.956194383870642;
                params(end).group = {'Q2'};
                params(end).label = 'beta Heterogeneity';
                
                % if ibw == 0.005
                %     params(end).betaH0 = -1e-3;
                % elseif ibw == 0.01
                %     params(end).betaH0 = -1e-3;
                % end
            end

            % random beta heterogeneity
            for ibw = [0.01]
                for bs = [1/50, 1/10]
                    name = sprintf('Beta5, pSwitch%g, pSpacing%g', bs, ibw);
                    params(end+1) = setup.Params(ifreq, name, IncomeProcess);
                    params(end) = set_shared_fields(params(end), income_params);
                    params(end).nbeta = 5;
                    params(end).betawidth = ibw;
                    params(end).prob_zswitch = bs;
                    params(end).dieprob = deathp;
                    params(end).group = {'Q2'};
                    params(end).label = 'beta Heterogeneity';

                    if bs == 1 / 50
                        params(end).betaH0 = -1e-2;
                    else
                        params(end).betaH0 = -1e-3;
                    end
                end
            end
        end

        % Different risk aversion coeffs
        for ira = [0.5, 2, 6]
            name = sprintf('CRRA = %g', ira);
            params(end+1) = setup.Params(4, name, quarterly_b_path);
            params(end) = set_shared_fields(params(end), quarterly_b_params);
            params(end).risk_aver = ira;
            if (ifreq==4 && ira==4) || ira==6
                params(end).betaL = 0.5;
            end

            if ira == 6
                params(end).betaH0 = -1e-3;
            end
            params(end).group = {'Q3'};
            params(end).label = 'CRRA';
        end
        
        % CRRA with IES heterogeneity
        name = 'CRRA w/IES betw exp(-1), exp(1)';
        params(end+1) = setup.Params(ifreq, name, IncomeProcess);
        params(end) = set_shared_fields(params(end), income_params);
        params(end).risk_aver = 1./ exp([-1 -0.5 0 0.5 1]);
        if params(end).freq == 4
            params(end).betaH0 =  - 5e-3;
        end
        params(end).group = {'Q3'};
        params(end).label = 'CRRA';
        params(end).other = {'exp(1), ..., exp(-1)', 'exp(-1), ..., exp(1)'};
        
        name = 'CRRA w/IES betw exp(-2), exp(2)';
        params(end+1) = setup.Params(ifreq, name, IncomeProcess);
        params(end) = set_shared_fields(params(end), income_params);
        params(end).risk_aver = 1./ exp([-2 -1 0 1 2]);
        if params(end).freq == 4
            params(end).betaH0 = -5e-3;
        end
        params(end).beta0 = 0.911905140057402;
        params(end).group = {'Q3'};
        params(end).label = 'CRRA';
        params(end).other = {'exp(2), ..., exp(-2)', 'exp(-2), ..., exp(2)'};

        % epstein-zin, quarterly
        ras = [0.5 8  1    1 8];
        ies = [1   1  0.25 2 2];
        for i = 1:5
            ra_i = ras(i);
            ies_i = ies(i);
            name = sprintf('EZ, ra%g, ies%g', ra_i, ies_i);
            params(end+1) = setup.Params(4, name, IncomeProcess);
            params(end) = set_shared_fields(params(end), income_params);
            params(end).risk_aver = ra_i;
            params(end).invies = 1 / ies_i;
            params(end).EpsteinZin = true;
            params(end).group = {'Q4'};
            params(end).label = 'EZ';

            if i <= 3
                params(end).betaH0 = - 3e-3;
            else
                params(end).betaH0 = - 1e-3;
            end
        end

        % EZ with IES heterogeneity
        name = 'EZ w/ IES betw exp(-1), exp(1)';
        params(end+1) = setup.Params(ifreq, name, IncomeProcess);
        params(end) = set_shared_fields(params(end), income_params);
        params(end).invies = 1 ./ exp([-1 -0.5 0 0.5 1]);
        params(end).EpsteinZin = true;
        if (ifreq == 4)
            params(end).betaH0 = - 3e-3;
        end
        params(end).group = {'Q4'};
        params(end).label = 'EZ';
        params(end).other = {1, 'exp(-1), ..., exp(1)'};
        
        name = 'EZ w/ IES betw exp(-2), exp(2)';
        params(end+1) = setup.Params(ifreq, name, IncomeProcess);
        params(end) = set_shared_fields(params(end), income_params);
        params(end).invies = 1 ./ exp([-2 -1 0 1 2]);
        params(end).EpsteinZin = true;
        if (ifreq == 4)
            params(end).betaH0 = - 3e-3;
        end
        params(end).group = {'Q4'};
        params(end).label = 'EZ';
        params(end).other = {1, 'exp(-2), ..., exp(2)'};

        % EZ with risk aversion heterogeneity
        name = 'EZ w/ RA betw exp(-2), exp(2)';
        params(end+1) = setup.Params(ifreq, name, IncomeProcess);
        params(end) = set_shared_fields(params(end), income_params);
        params(end).invies = 1;
        params(end).risk_aver = exp([-2 -1 0 1 2]);
        params(end).EpsteinZin = 1;
        params(end).betaH0 = - 3e-3;
        params(end).betaL = 0.96;
        params(end).beta0 = 0.99^4;
        params(end).group = {'Q4'};
        params(end).label = 'EZ';
        params(end).other = {'exp(-2), ..., exp(2)', 1};
    end

    % temptation
    for tempt = [0.01 0.05 0.07]
        name = sprintf('Temptation = %g', tempt);
        params(end+1) = setup.Params(4, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end).temptation = tempt;
        if tempt == 0.07
            params(end).betaH0 = 6e-4;
        elseif  tempt == 0.05
            params(end).betaH0 = -2e-5;
        end
        params(end).group = {'Q5'};
        params(end).label = 'Temptation';
    end

    name = 'Temptation, uniform in {0, 0.05, 0.1}';
    params(end+1) = setup.Params(4, name, quarterly_b_path);
    params(end) = set_shared_fields(params(end), quarterly_b_params);
    params(end).temptation = [0 0.05 0.1];
    params(end).beta0 = 0.998283 ^ 4;
    params(end).betaH0 = 5e-2;
    params(end).group = {'Q5'};
    params(end).label = 'Temptation';
    params(end).other = {'{0, 0.05, 0.1}'};

    %----------------------------------------------------------------------
    % PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    name = 'Annual, no yT';
    params(end+1) = setup.Params(1, name, '');
    params(end) = set_shared_fields(params(end), annual_params);
    params(end).beta0 = 0.99;
    params(end).nyT = 1;
    params(end).sd_logyT = 0;
    params(end).lambdaT = 0;
    params(end).group = {'Q7'};
    params(end).label = 'Annual (i)';
    params(end).other = {'no trans shocks'};

%     % ii
%     params(end+1) = setup.Params(1,'A a(ii) MeasError','');
%     params(end).sd_logyT = sqrt(0.02);
%     
%     % iii
%     params(end+1) = setup.Params(1,'A a(iii) NoTranReEst','');
%     params(end).rho_logyP = 0.8592;
%     params(end).sd_logyP = sqrt(0.132);
%     params(end).nyT = 1;
%     params(end).sd_logyT = 0;

    % iv
    name = 'Annual, Carrol';
    params(end+1) = setup.Params(1, name, '');
    params(end) = set_shared_fields(params(end), annual_params);
    params(end).rho_logyP = 0.999;
    params(end).sd_logyP = sqrt(0.015);
    params(end).sd_logyT = sqrt(0.01);
    params(end).group = {'Q7'};
    params(end).label = 'Annual (iv)';
    params(end).other = {'Carrol process'};
    
%     % v
%     params(end+1) = setup.Params(1,'A a(v) HighPersNotReEst','');
%     params(end).rho_logyP = 0.99;
%     
%     % vi
%     params(end+1) = setup.Params(1,'A a(vi) LowPersNotReEst','');
%     params(end).rho_logyP = 0.9;
% 
%     % vii
%     params(end+1) = setup.Params(1,'A a(vii) HighPersReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).sd_logyT = sqrt(0.0667);
    
    % viii
    name = 'Annual, high persistence';
    params(end+1) = setup.Params(1, name, '');
    params(end) = set_shared_fields(params(end), shared_params);
    params(end).rho_logyP = 0.995;
    params(end).sd_logyP = sqrt(0.0043);
    params(end).sd_logyT = sqrt(0.0688);
    params(end).lambdaT = 1;
    params(end).group = {'Q7'};
    params(end).label = 'Annual (viii)';
    params(end).other = {'high persistence'};
%     
%     % ix
%     params(end+1) = Params(1,'A a(ix) HighPersNoTransReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).nyT = 1;
%     params(end).sd_logyT = sqrt(0);
    
    % x
    name = 'Annual, high nyF = 5';
    params(end+1) = setup.Params(1, name, '');
    params(end) = set_shared_fields(params(end), shared_params);
    params(end).rho_logyP = 0.9158;
    params(end).sd_logyP = sqrt(0.0445);
    params(end).sd_logyT = sqrt(0.0479);
    params(end).sd_logyF = sqrt(0.1801);
    params(end).lambdaT = 1;
    params(end).nyF = 5;
    params(end).group = {'Q7'};
    params(end).label = 'Annual (x)';
    params(end).other = {'FE heterogeneity'};

%     % xi
%     params(end+1) = setup.Params(1,'A a(xi) MatchSSA','');
%     params(end).rho_logyP = 0.9468;
%     params(end).sd_logyP = sqrt(0.0641);
%     params(end).sd_logyT = sqrt(0.0479);
%     params(end).lambdaT  = 0.0821;
%     
%     % xii
%     params(end+1) = setup.Params(1,'A a(xii) WithSCF m0','');
%     params(end).rho_logyP = 0.9787;
%     params(end).sd_logyP = sqrt(0.0400);
%     params(end).sd_logyT = sqrt(0.0508);
%     
%     % xiv
%     params(end+1) = setup.Params(1,'A a(xiv) MassPointTrans','');
%     params(end).rho_logyP = sqrt(0.9516);
%     params(end).sd_logyP = sqrt(0.0434);
%     params(end).sd_logyT = sqrt(0.6431);
%     params(end).lambdaT = 0.0760;
    
    %----------------------------------------------------------------------
    % PART 3b, QUARTERLY MODEL
    %----------------------------------------------------------------------
    
    
    % i quarterly_a
    name = 'quarterly_a';
    params(end+1) = setup.Params(4, name , '');
    params(end) = set_shared_fields(params(end), quarterly_a_params);
    params(end).beta0 = 0.984363510593659;
    params(end).group = {'Q7'};
    params(end).label = 'Quart (i)';
    params(end).other = {'quart_a'};
    
    % ii
    name = 'KMP';
    params(end+1) = setup.Params(4, name, '');
    params(end) = set_shared_fields(params(end), shared_params);
    params(end).rho_logyP = 0.9879;
    params(end).sd_logyP = sqrt(0.0109);
    params(end).sd_logyT = sqrt(0.0494);
    params(end).group = {'Q7'};
    params(end).label = 'Quart (ii)';
    params(end).other = {'KMP'};

    % % KMP with tax and transfer - Mitman inc process
    % params(end+1) = setup.Params(4, 'Q KMP (Mitman income) w/tax and transfer, no discount het', 'input/income_mitman.mat');
    % params(end).sd_logyT = sqrt(0.0522);
    % params(end).labtaxlow = 0.25;
    % params(end).lumptransfer = 0.0363;
    % params(end).target_value = 4.9;
    % params(end).r = 0;
    
    
    % % KMP with tax and transfer - our inc process
    % params(end+1) = setup.Params(4, 'Q KMP (our income) w/tax and transfer, no discount het', '');
    % params(end).rho_logyP = 0.9879;
    % params(end).sd_logyP = sqrt(0.0109);
    % params(end).sd_logyT = sqrt(0.0494);
    % params(end).lambdaT = 1;
    % params(end).labtaxlow = 0.25;
    % params(end).lumptransfer = 0.0363;
    % params(end).target_value = 4.9;
    % params(end).r = 0;
    
    % % IMP with tax and transfer, and discount factor heterogeneity- Mitman inc process
    % params(end+1) = setup.Params(4, 'Q KMP (Mitman income) w/tax and transfer, beta width 0.01', 'input/mitman.mat');
    % params(end).sd_logyT = sqrt(0.0522);
    % params(end).lambdaT = 1;
    % params(end).labtaxlow = 0.25;
    % params(end).lumptransfer = 0.0363;
    % params(end).target_value = 4.9;
    % params(end).r = 0;
    % params(end).nbeta = 2;
    % params(end).betawidth = 0.01;
    % params(end).beta_dist = [0.2, 0.8];
    % params(end).beta0 = 0.9;
    % params(end).betaH0 = -1e-3;

    %  % IMP with tax and transfer, and discount factor heterogeneity- Mitman inc process
    % params(end+1) = setup.Params(4, 'Q KMP (Mitman income) w/tax and transfer, beta width 0.1', 'input/mitman.mat');
    % params(end).sd_logyT = sqrt(0.0522);
    % params(end).lambdaT = 1;
    % params(end).labtaxlow = 0.25;
    % params(end).lumptransfer = 0.0363;
    % params(end).target_value = 4.9;
    % params(end).r = 0;
    % params(end).nbeta = 2;
    % params(end).betawidth = 0.1;
    % params(end).beta_dist = [0.2, 0.8];
    % params(end).beta0 = 0.9;
    % params(end).betaH0 = -1e-3;

    % % IMP with tax and transfer, and discount factor heterogeneity- Mitman inc process
    % params(end+1) = setup.Params(4, 'Q KMP (Mitman income) w/tax and transfer, beta 0.9929, 0.9994', 'input/mitman.mat');
    % params(end).sd_logyT = sqrt(0.0522);
    % params(end).lambdaT = 1;
    % params(end).labtaxlow = 0.25;
    % params(end).lumptransfer = 0.0363;
    % params(end).r = 0;
    % params(end).nbeta = 2;
    % params(end).beta_dist = [0.2, 0.8];
    % params(end).beta_grid_forced = [0.9929; 0.9994];
    
    
    % iii quarterly_c
    name = 'quarterly_c';
    params(end+1) = setup.Params(4, name, quarterly_c_path);
    params(end) = set_shared_fields(params(end), quarterly_c_params);
    params(end).beta0 = 0.984363510593659;
    params(end).group = {'Q7'};
    params(end).label = 'Quart (iii)';
    params(end).other = {'quart_c'};

%     % iv
%     params(end+1) = setup.Params(4,'Q b(iv) PersEveryPeriod','');
%     params(end).rho_logyP = 0.9884;
%     params(end).sd_logyP = sqrt(0.0105);
%     params(end).sd_logyT = sqrt(1.5298);
%     params(end).lambdaT = 0.0813;

    %----------------------------------------------------------------------
    % PART 4, Exotic Preferences
    %----------------------------------------------------------------------
    for ifreq = [4]
        if ifreq == 1
            lfreq = 'A';
            IncomeProcess = '';
            income_params = annual_params;
        else
            lfreq = 'Q';
            IncomeProcess = quarterly_b_path;
            income_params = quarterly_b_params;
        end
        
        
    end
        
    
        
%         % epstein-zin: vary risk_aver
%         for ra = [0.5 0.75 1.5 2 4 8]
%             params(end+1) = setup.Params(ifreq,['EZ ra' num2str(ra) ' ies1'],IncomeProcess);
%             params(end).risk_aver = ra;
%             params(end).invies = 1;
%             params(end).EpsteinZin = 1;
%         end
%         
%         % epstein-zin: vary invies
%         for ies = [1/4 1/2 3/4 1.5 2 5]
%             params(end+1) = setup.Params(ifreq,['EZ ra1 ies' num2str(ies)],IncomeProcess);
%             params(end).risk_aver = 1;
%             params(end).invies = 1/ies;
%             params(end).EpsteinZin = 1;
%         end

    %----------------------------------------------------------------------
    % OTHER
    %----------------------------------------------------------------------
    % params(end+1) = setup.Params(4, 'quarterly_b_nyT101', quarterly_b_path);
    % params(end) = set_shared_fields(params(end), quarterly_b_params);
    % params(end).nyT = 101;
    % params(end).beta0 = 0.984363510593659;

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS, DO NOT CHANGE
    %----------------------------------------------------------------------
    params.set_index();

    % get list of all names
    all_names = cell2table({params.name}');
    
    % select by number if there is one, otherwise select by names,
    % otherwise use all
    if numel(runopts.number) == 1
        params = params(runopts.number);
    elseif numel(runopts.number) > 1
        error('runopts.number must have 1 or zero elements')
    else
        params = setup.Params.select_by_names(params, runopts.names_to_run);
    end

    params.set_run_parameters(runopts);
    params.make_adjustments();

    %----------------------------------------------------------------------
    % ATTACH CALIBRATOR
    %----------------------------------------------------------------------
    if params.calibrate
        heterogeneity = setup.Prefs_R_Heterogeneity(params);

        if (params.nbeta > 1) && isequal(heterogeneity.ztrans, eye(params.nbeta))
            new_betaH = params.betaH - max(heterogeneity.betagrid0);
            params.set("betaH", new_betaH, true);
        end

        calibrator = aux.mean_wealth_calibrator(params);
        calibrator.set_handle(params);
        params.set("calibrator", calibrator, true);
    end
end
