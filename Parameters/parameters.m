function params = parameters(runopts,selection,QIncome)
    
    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    
    % Annual
    params(1) = MPCParams(1,'baseline_A','');
    
    % Quarterly
    params(end+1) = MPCParams(4,'baseline_Q',QIncome);
    
    %----------------------------------------------------------------------
    % PART 2, DIFFERENT ASSUMPTIONS
    %----------------------------------------------------------------------

    for ifreq = [1 4]
        if ifreq == 1
            lfreq = 'A';
            IncomeProcess = '';
        else
            lfreq = 'Q';
            IncomeProcess = QIncome;
        end
        % different mean wealth targets
        for mw = [0.25, 0.5, 1]
            name = [lfreq ' AYtarget' num2str(mw) ];
            params(end+1) = MPCParams(ifreq,name,IncomeProcess);
            params(end).targetAY = mw;
            if ifreq == 4
                params(end).betaL = 0.5;
            end
        end

        % different interest rates
        for ii = [0, 5]
            name = [lfreq ' IntRate' num2str(ii)];
            params(end+1) = MPCParams(ifreq,name,IncomeProcess);
            params(end).r = ii/100;
        end

%         % different tax rates
%         for itax = [0.05, 0.1, 0.15, 0.25]
%             name = [lfreq ' LabTax' num2str(itax)];
%             params(end+1) = MPCParams(ifreq,name,'');
%             params(end).labtaxlow = itax;
%         end

        % no death
        name = [lfreq ' NoDeath'];
        params(end+1) = MPCParams(ifreq,name,IncomeProcess);
        params(end).dieprob = 0;

        % no bequests
        name = [lfreq ' NoBequests'];
        params(end+1) = MPCParams(ifreq,name,IncomeProcess);
        params(end).Bequests = 0;

        % perfect annuities
        name = [lfreq ' Annuities'];
        params(end+1) = MPCParams(ifreq,name,IncomeProcess);
        params(end).annuities_on();

%         % bequest curvature
%         for bcurv = [0.1 0.5 1 2 5]
%             name = [lfreq ' BeqWt0.02 BeqLux0.01 BeqCurv' num2str(bcurv)];
%             params(end+1) = MPCParams(ifreq,name,IncomeProcess);
%             params(end).bequest_weight = 0.02;
%             params(end).bequest_luxury = 0.01;
%             params(end).bequest_curv   = bcurv;
%         end

       
        for deathp = 1/50 %[0 1/50]
            if deathp == 0
                deathind = ' NoDeath';
            else
                deathind = ' Death';
            end
             % fixed beta heterogeneity
            for ibw = [0.001, 0.005, 0.01]
                name = [lfreq ' FixedBetaHet5 Width' num2str(ibw) deathind];
                params(end+1) = MPCParams(ifreq,name,IncomeProcess);
                params(end).nb = 5;
                params(end).betawidth = ibw;
                params(end).betaswitch = 0;
                params(end).dieprob = deathp;
            end

            % random beta heterogeneity
            for ibw = [0.01]
                for bs = [1/50, 1/10]
                    name = [lfreq ' RandomBetaHet5 Width' num2str(ibw) ' SwitchProb' num2str(bs) deathind];
                    params(end+1) = MPCParams(ifreq,name,IncomeProcess);
                    params(end).nb = 5;
                    params(end).betawidth = ibw;
                    params(end).betaswitch = bs;
                    params(end).dieprob = deathp;
                end
            end
        end
        
        % CRRA with IES heterogeneity
        params(end+1) = MPCParams(ifreq,[lfreq ' CRRA with IES heterogeneity'],IncomeProcess);
        params(end).risk_aver = exp([-1 -0.5 0 0.5 1]);

        % EZ with IES heterogeneity
        params(end+1) = MPCParams(ifreq,[lfreq ' EZ with IES heterogeneity'],IncomeProcess);
        params(end).invies = 1 ./ exp([-1 -0.5 0 0.5 1]);
        params(end).EpsteinZin = 1;
    end

    %----------------------------------------------------------------------
    % PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    params(end+1) = MPCParams(1,'A a(i) NoTransShocks','');
    params(end).nyT = 1;
    params(end).sd_logyT = 0;
    
%     % ii
%     params(end+1) = MPCParams(1,'A a(ii) MeasError','');
%     params(end).sd_logyT = sqrt(0.02);
%     
%     % iii
%     params(end+1) = MPCParams(1,'A a(iii) NoTranReEst','');
%     params(end).rho_logyP = 0.8592;
%     params(end).sd_logyP = sqrt(0.132);
%     params(end).nyT = 1;
%     params(end).sd_logyT = 0;

    % iv
    params(end+1) = MPCParams(1,'A a(iv) HighPersistCarrol','');
    params(end).rho_logyP = 0.999;
    params(end).sd_logyP = sqrt(0.015);
    params(end).sd_logyT = sqrt(0.01);
    
%     % v
%     params(end+1) = MPCParams(1,'A a(v) HighPersNotReEst','');
%     params(end).rho_logyP = 0.99;
%     
%     % vi
%     params(end+1) = MPCParams(1,'A a(vi) LowPersNotReEst','');
%     params(end).rho_logyP = 0.9;
% 
%     % vii
%     params(end+1) = MPCParams(1,'A a(vii) HighPersReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).sd_logyT = sqrt(0.0667);
    
    % viii
    params(end+1) = MPCParams(1,'A a(viii) EvenHigherPersReEst','');
    params(end).rho_logyP = 0.995;
    params(end).sd_logyP = sqrt(0.0043);
    params(end).sd_logyT = sqrt(0.0688);
%     
%     % ix
%     params(end+1) = MPCParams(1,'A a(ix) HighPersNoTransReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).nyT = 1;
%     params(end).sd_logyT = sqrt(0);
    
    % x
    params(end+1) = MPCParams(1,'A WithFE nyF 5','');
    params(end).rho_logyP = 0.9158;
    params(end).sd_logyP = sqrt(0.0445);
    params(end).sd_logyT = sqrt(0.0479);
    params(end).sd_logyF = sqrt(0.1801);
    params(end).nyF = 5;

%     % xi
%     params(end+1) = MPCParams(1,'A a(xi) MatchSSA','');
%     params(end).rho_logyP = 0.9468;
%     params(end).sd_logyP = sqrt(0.0641);
%     params(end).sd_logyT = sqrt(0.0479);
%     params(end).lambdaT  = 0.0821;
%     
%     % xii
%     params(end+1) = MPCParams(1,'A a(xii) WithSCF m0','');
%     params(end).rho_logyP = 0.9787;
%     params(end).sd_logyP = sqrt(0.0400);
%     params(end).sd_logyT = sqrt(0.0508);
%     
%     % xiv
%     params(end+1) = MPCParams(1,'A a(xiv) MassPointTrans','');
%     params(end).rho_logyP = sqrt(0.9516);
%     params(end).sd_logyP = sqrt(0.0434);
%     params(end).sd_logyT = sqrt(0.6431);
%     params(end).lambdaT = 0.0760;
    
    %----------------------------------------------------------------------
    % PART 3b, QUARTERLY MODEL
    %----------------------------------------------------------------------
    
    % different risk aversion coeffs
    for ira = [0.5, 2, 6]
        name = ['Q RiskAver' num2str(ira)];
        params(end+1) = MPCParams(4,name,QIncome);
        params(end).risk_aver = ira;
        if (ifreq==4 && ira==4) || ira==6
            params(end).betaL = 0.5;
        end
    end
    
    % i quarterly_a
    params(end+1) = MPCParams(4,'Q b(i) quarterly_a','');
    
    % ii
    params(end+1) = MPCParams(4,'Q b(ii) KMPTransf','');
    params(end).rho_logyP = 0.9879;
    params(end).sd_logyP = sqrt(0.0109);
    params(end).sd_logyT = sqrt(0.0494);
    
    % iii quarterly_c
    params(end+1) = MPCParams(4,'Q b(iii) quarterly_c','IncomeGrids/quarterly_c.mat');
    
%     % iv
%     params(end+1) = MPCParams(4,'Q b(iv) PersEveryPeriod','');
%     params(end).rho_logyP = 0.9884;
%     params(end).sd_logyP = sqrt(0.0105);
%     params(end).sd_logyT = sqrt(1.5298);
%     params(end).lambdaT = 0.0813;

    %----------------------------------------------------------------------
    % PART 4, Exotic Preferences
    %----------------------------------------------------------------------
    for ifreq = [1 4]
        if ifreq == 1
            lfreq = 'A';
            IncomeProcess = '';
        else
            lfreq = 'Q';
            IncomeProcess = QIncome;
        end
        
        % temptation
        for itempt = [0.005 0.01 0.05]
            params(end+1) = MPCParams(ifreq,[lfreq ' Temptation' num2str(itempt)],IncomeProcess);
            params(end).temptation = itempt;
            if itempt == 0.05 && ifreq == 4
                params(end).set_betaH_distance(-1e-5);
            end    
        end
    end
        
    % epstein-zin, quarterly
    ras = [0.5 8  1    1 8];
    ies = [1   1  0.25 2 2];
    for i = 1:5
        params(end+1) = MPCParams(4,['Q EZ ra' num2str(ras(i)) ' ies' num2str(ies(i))],QIncome);
        params(end).risk_aver = ras(i);
        params(end).invies = 1 / ies(i);
        params(end).EpsteinZin = 1;
    end
        
%         % epstein-zin: vary risk_aver
%         for ra = [0.5 0.75 1.5 2 4 8]
%             params(end+1) = MPCParams(ifreq,['EZ ra' num2str(ra) ' ies1'],IncomeProcess);
%             params(end).risk_aver = ra;
%             params(end).invies = 1;
%             params(end).EpsteinZin = 1;
%         end
%         
%         % epstein-zin: vary invies
%         for ies = [1/4 1/2 3/4 1.5 2 5]
%             params(end+1) = MPCParams(ifreq,['EZ ra1 ies' num2str(ies)],IncomeProcess);
%             params(end).risk_aver = 1;
%             params(end).invies = 1/ies;
%             params(end).EpsteinZin = 1;
%         end

    %----------------------------------------------------------------------
    % ADJUST TO QUARTERLY VALUES
    %----------------------------------------------------------------------

    params = MPCParams.adjust_if_quarterly(params);
    
    %----------------------------------------------------------------------
    % SET BETA UPPER BOUND FOR SPECIAL CASES
    %----------------------------------------------------------------------
    % This section adjusts the beta upper bound for specific cases
    % to ensure convergence

    % CRRA heterogeneity
    params.set_betaH_distance(-2e-2,'A CRRA with IES heterogeneity',1);
    params.set_betaH_distance(-5e-3,'Q CRRA with IES heterogeneity',4);

    % temptation
    params.set_betaH_distance(-1e-5,'Q Temptation0.05',4);
    params.set_betaH_distance(-1e-5,'A Temptation0.05',1);
    
%     % Epstein-Zin
    params.set_betaH_distance(-3e-2,'A EZ with IES heterogeneity',1);
    params.set_betaH_distance(-8e-3,'Q EZ with IES heterogeneity',4);
    
%     params.set_betaH_distance(-3e-2,'A EZ ra0.5 ies1',1);
    params.set_betaH_distance(-8e-3,'Q EZ ra0.5 ies1',4);
%     params.set_betaH_distance(-3e-2,'A EZ ra8 ies1',1);
    params.set_betaH_distance(-8e-3,'Q EZ ra8 ies1',4);
%     params.set_betaH_distance(-3e-2,'A EZ ra1 ies0.25',1);
    params.set_betaH_distance(-8e-3,'Q EZ ra1 ies0.25',4);
%     params.set_betaH_distance(-2.5e-2,'A EZ ra1 ies2',1);
    params.set_betaH_distance(-6.5e-3,'Q EZ ra1 ies2',4);
%     params.set_betaH_distance(-2.5e-2,'A EZ ra8 ies2',1);
    params.set_betaH_distance(-6.5e-3,'Q EZ ra8 ies2',4);

%     % varying risk_aver
%     EZ = find([params.EpsteinZin]==1 & [params.freq]==1);
%     params.set_betaH_distance(-3e-2,name,1);
% 
%     EZ = find([params.EpsteinZin]==1 & [params.freq]==4);
%     params.set_betaH_distance(-6e-3,name,4);
    
    
%     % --------- annual, varying invies -------------
%     freq = 1;
%     
%     change_betaH = ['EZ ra1 ies' num2str(1/4)];
%     params.set_betaH_distance(-3e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(1/2)];
%     params.set_betaH_distance(-3e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(3/4)];
%     params.set_betaH_distance(-3e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(1.5)];
%     params.set_betaH_distance(-2.1e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(2)];
%     params.set_betaH_distance(-2.1e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(5)];
%     params.set_betaH_distance(-2.1e-2,change_betaH,freq);
%     
%     % --------- quarterly, varying invies -------------
%     freq = 4;
%     
%     change_betaH = ['EZ ra1 ies' num2str(1/4)];
%     params.set_betaH_distance(-1.5e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(1/2)];
%     params.set_betaH_distance(-1.5e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(3/4)];
%     params.set_betaH_distance(-1e-2,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(1.5)];
%     params.set_betaH_distance(-5.5e-3,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(2)];
%     params.set_betaH_distance(-5.5e-3,change_betaH,freq);
%     
%     change_betaH = ['EZ ra1 ies' num2str(5)];
%     params.set_betaH_distance(-5.4e-3,change_betaH,freq);
    
    
    %----------------------------------------------------------------------
    % SET BETA UPPER BOUND FOR BETA HETEROGENEITY CASES
    %----------------------------------------------------------------------
    
    change_betaH = [' RandomBetaHet5 Width' num2str(0.01)...
                        ' SwitchProb' num2str(1/10) ' Death'];
    params.set_betaH_distance(1e-2,['A' change_betaH],1);
    params.set_betaH_distance(1.3e-2,['Q' change_betaH],4);
    
    change_betaH = [' RandomBetaHet5 Width' num2str(0.01)...
                        ' SwitchProb' num2str(1/50) ' Death'];
    params.set_betaH_distance(1e-2,['A' change_betaH],1);
    params.set_betaH_distance(5e-3,['Q' change_betaH],4);
    
%     % --------- annual, random ----------
%     freq = 1;
% 
%     % no death
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 NoDeath';
%     params.set_betaH_distance(-3e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 NoDeath';
%     params.set_betaH_distance(-3e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 NoDeath';
%     params.set_betaH_distance(-3e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 NoDeath';
%     params.set_betaH_distance(-1e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 NoDeath';
%     params.set_betaH_distance(-8e-4,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 NoDeath';
%     params.set_betaH_distance(5e-3,change_betaH,freq);
% 
%     % death
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 Death';
%     params.set_betaH_distance(-1.3e-2,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 Death';
%     params.set_betaH_distance(-9e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 Death';
%     params.set_betaH_distance(-7e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 Death';
%     params.set_betaH_distance(-4e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 Death';
%     params.set_betaH_distance(-1e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 Death';
%     params.set_betaH_distance(3e-3,change_betaH,freq);
% 
%     % --------- quarterly, random -----------
%     freq = 4;
% 
%     % no death
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 NoDeath';
%     params.set_betaH_distance(-1e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 NoDeath';
%     params.set_betaH_distance(-1e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 NoDeath';
%     params.set_betaH_distance(3e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 NoDeath';
%     params.set_betaH_distance(5e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 NoDeath';
%     params.set_betaH_distance(4e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 NoDeath';
%     params.set_betaH_distance(1.5e-2,change_betaH,freq);
% 
%     % death
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 Death';
%     params.set_betaH_distance(-1e-4,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 Death';
%     params.set_betaH_distance(-1e-4,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 Death';
%     params.set_betaH_distance(2e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 Death';
%     params.set_betaH_distance(5.1e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 Death';
%     params.set_betaH_distance(4.7e-3,change_betaH,freq);
% 
%     change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 Death';
%     params.set_betaH_distance(1.2e-2,change_betaH,freq);
    
    %----------------------------------------------------------------------
    % SET BETA UPPER BOUND FOR OTHER CASES
    %----------------------------------------------------------------------
    
%     change_betaH = '2 BeqWt0.02 BeqLux0.01 BeqCurv0.1';
%     params.set_betaH_distance(-5e-3,['A' change_betaH],1);
%     change_betaH = ' Temptation0.05';
%     params.set_betaH_distance(-1e-4,['Q' change_betaH],4);S

    change_betaH = ' RiskAver6';
    params.set_betaH_distance(-1e-2,['Q' change_betaH],4);

    change_betaH = ' Annuities';
    params.set_betaH_distance(-5e-3,['Q' change_betaH],4);

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------

    params.set_run_parameters(runopts);

    % creates ordered 'index' field
    params.set_index();
    
    % select by number if there is one, otherwise select by names,
    % otherwise use all
    if numel(selection.number) == 1
        params = MPCParams.select_by_number(params,selection.number);
    elseif numel(selection.number) > 1
        error('selection.number must have 1 or zero elements')
    else
        params = MPCParams.select_by_names(params,selection.names_to_run);
        params.set_index(); % index within .mat file
    end
   
    

end
