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

    for ifreq = [4]
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
        params(end).betaH0 = params(end).betaH0 - 5e-3;

%         % bequest curvature
%         for bcurv = [0.1 0.5 1 2 5]
%             name = [lfreq ' BeqWt0.02 BeqLux0.01 BeqCurv' num2str(bcurv)];
%             params(end+1) = MPCParams(ifreq,name,IncomeProcess);
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
                    if (bs == 1/50) && (deathp == 1/50)
                        params(end).betaH0 = params(end).betaH0 + 1.3e-2;
                    elseif (deathp == 1/50)
                        params(end).betaH0 = params(end).betaH0 + 5e-3;
                    end
                end
            end
        end

        
        % CRRA with IES heterogeneity
        params(end+1) = MPCParams(ifreq,[lfreq ' CRRA w/IES betw exp(-1), exp(1)'],IncomeProcess);
        params(end).risk_aver = 1./ exp([-1 -0.5 0 0.5 1]);
        if params(end).freq == 4
            params(end).betaH0 =  params(end).betaH0 - 5e-3;
        end
        
        params(end+1) = MPCParams(ifreq,[lfreq ' CRRA w/IES betw exp(-2), exp(2)'],IncomeProcess);
        params(end).risk_aver = 1./ exp([-2 -1 0 1 2]);
        if params(end).freq == 4
            params(end).betaH0 =  params(end).betaH0 - 5e-3;
        end

        % EZ with IES heterogeneity
        params(end+1) = MPCParams(ifreq,[lfreq ' EZ w/IES betw exp(-1), exp(1)'],IncomeProcess);
        params(end).invies = 1 ./ exp([-1 -0.5 0 0.5 1]);
        params(end).EpsteinZin = 1;
        if (ifreq == 4)
            params(end).betaH0 = params(end).betaH0 - 8e-3;
        end
        
        params(end+1) = MPCParams(ifreq,[lfreq ' EZ w/IES betw exp(-2), exp(2)'],IncomeProcess);
        params(end).invies = 1 ./ exp([-2 -1 0 1 2]);
        params(end).EpsteinZin = 1;
        if (ifreq == 4)
            params(end).betaH0 = params(end).betaH0 - 8e-3;
        end
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

        if ira == 6
            params(end).betaH0 = params(end).betaH0 - 1e-2;
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
    for ifreq = [4]
        if ifreq == 1
            lfreq = 'A';
            IncomeProcess = '';
        else
            lfreq = 'Q';
            IncomeProcess = QIncome;
        end
        
        % temptation
        for itempt = [0.01 0.05 0.07]
            params(end+1) = MPCParams(ifreq,[lfreq ' Temptation' num2str(itempt)],IncomeProcess);
            params(end).temptation = itempt;
            if (ifreq==4) && (itempt==0.07)
                params(end).betaH0 = params(end).betaH0 + 3.2e-4;
            elseif (ifreq==4) && (itempt==0.05)
                params(end).betaH0 = params(end).betaH0 - 2e-5;
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
        if i <= 3
            params(end).betaH0 = params(end).betaH0 - 8e-3;
        else
            params(end).betaH0 = params(end).betaH0 - 6.5e-3;
        end
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
