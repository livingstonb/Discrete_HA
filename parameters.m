function params = parameters(runopts,selection)
    
    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    
    % Annual
    params(1) = MPCParams(1,'baseline_A');
    
    % Quarterly
    params(end+1) = MPCParams(4,'baseline_Q');
    
    %----------------------------------------------------------------------
    % PART 2, DIFFERENT ASSUMPTIONS
    %----------------------------------------------------------------------
    for ifreq = [1 4]
        % different mean wealth targets
        for mw = [0.25, 0.5, 1, 2.5, 5]
            name = ['2 AYtarget' num2str(mw) ];
            params(end+1) = MPCParams(ifreq,name);
            params(end).targetAY = mw;
            if ifreq == 4
                params(end).betaL = 0.5;
            end
        end

        % different interest rates
        for ii = [0, 1, 3, 5]
            name = ['2 IntRate' num2str(ii)];
            params(end+1) = MPCParams(ifreq,name);
            params(end).r = ii/100;
        end

        % different risk aversion coeffs
        for ira = [0.5, 1.5, 2, 4, 6]
            name = ['2 RiskAver' num2str(ira)];
            params(end+1) = MPCParams(ifreq,name);
            params(end).risk_aver = ira;
            if (ifreq==4 && ira==4) || ira==6
                params(end).betaL = 0.5;
            end
        end

        % different tax rates
        for itax = [0.05, 0.1, 0.15, 0.25]
            name = ['2 LabTax' num2str(itax)];
            params(end+1) = MPCParams(ifreq,name);
            params(end).labtaxlow = itax;
        end

        % no death
        params(end+1) = MPCParams(ifreq,'2 NoDeath');
        params(end).dieprob = 0;

        % no bequests
        params(end+1) = MPCParams(ifreq,'2 NoBequests');
        params(end).Bequests = 0;

        % perfect annuities
        params(end+1) = MPCParams(ifreq,'2 Annuities');
        params(end).annuities_on();

        % bequest curvature
        for bcurv = [0.1 0.5 1 2 5]
            name = ['2 BeqWt0.02 BeqLux0.01 BeqCurv' num2str(bcurv)];
            params(end+1) = MPCParams(ifreq,name);
            params(end).bequest_weight = 0.02;
            params(end).bequest_luxury = 0.01;
            params(end).bequest_curv   = bcurv;
        end

        % fixed beta heterogeneity
        for deathp = [0 1/50]
            if deathp == 0
                deathind = ' NoDeath';
            else
                deathind = ' Death';
            end

            for ibw = [0.001, 0.005, 0.01]
                name = ['2 FixedBetaHet5 Width' num2str(ibw) deathind];
                params(end+1) = MPCParams(ifreq,name);
                params(end).nb = 5;
                params(end).nxlong = 400;
                params(end).nx = 100;
                params(end).betawidth = ibw;
                params(end).betaswitch = 0;
                params(end).dieprob = deathp;
            end

            % random beta heterogeneity
            for ibw = [0.001, 0.005, 0.01]
                for bs = [1/50, 1/10]
                    name = ['2 RandomBetaHet5 Width' num2str(ibw) ' SwitchProb' num2str(bs) deathind];
                    params(end+1) = MPCParams(ifreq,name);
                    params(end).nb = 5;
                    params(end).nxlong = 400;
                    params(end).nx = 100;
                    params(end).betawidth = ibw;
                    params(end).betaswitch = bs;
                    params(end).dieprob = deathp;
                end
            end
        end
    end

    %----------------------------------------------------------------------
    % PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    params(end+1) = MPCParams(1,'3a(i) NoTransShocks');
    params(end).nyT = 1;
    params(end).sd_logyT = 0;
    
    % ii
    params(end+1) = MPCParams(1,'3a(ii) MeasError');
    params(end).sd_logyT = sqrt(0.02);
    
    % iii
    params(end+1) = MPCParams(1,'3a(iii) NoTranReEst');
    params(end).rho_logyP = 0.8592;
    params(end).sd_logyP = sqrt(0.132);
    params(end).nyT = 1;
    params(end).sd_logyT = 0;

    % iv
%     params(end) = baseline;
%     params(end).name = '3a(iv) HighPersistCarrol';
%     params(end).rho_logyP = 0.9995;
%     params(end).sd_logyP = sqrt(0.015);
%     params(end).sd_logyT = sqrt(0.01);
    
    % v
    params(end+1) = MPCParams(1,'3a(v) HighPersNotReEst');
    params(end).rho_logyP = 0.99;
    
    % vi
    params(end+1) = MPCParams(1,'3a(vi) LowPersNotReEst');
    params(end).rho_logyP = 0.9;

    % vii
    params(end+1) = MPCParams(1,'3a(vii) HighPersReEst');
    params(end).rho_logyP = 0.99;
    params(end).sd_logyP = sqrt(0.0088);
    params(end).sd_logyT = sqrt(0.0667);
    
    % viii
    params(end+1) = MPCParams(1,'3a(viii) EvenHigherPersReEst');
    params(end).rho_logyP = 0.995;
    params(end).sd_logyP = sqrt(0.0043);
    params(end).sd_logyT = sqrt(0.0688);
    
    % ix
    params(end+1) = MPCParams(1,'3a(ix) HighPersNoTransReEst');
    params(end).rho_logyP = 0.99;
    params(end).sd_logyP = sqrt(0.0088);
    params(end).nyT = 1;
    params(end).sd_logyT = sqrt(0);
    
    % x
    params(end+1) = MPCParams(1,'WithFE nyF 3');
    params(end).rho_logyP = 0.9158;
    params(end).sd_logyP = sqrt(0.0445);
    params(end).sd_logyT = sqrt(0.0479);
    params(end).sd_logyF = sqrt(0.1801);
    params(end).nyF = 3;

    % xi
    params(end+1) = MPCParams(1,'3a(xi) MatchSSA');
    params(end).rho_logyP = 0.9468;
    params(end).sd_logyP = sqrt(0.0641);
    params(end).sd_logyT = sqrt(0.0479);
    params(end).lambdaT  = 0.0821;
    
    % xii
    params(end+1) = MPCParams(1,'3a(xii) WithSCF m0');
    params(end).rho_logyP = 0.9787;
    params(end).sd_logyP = sqrt(0.0400);
    params(end).sd_logyT = sqrt(0.0508);
    
    % xiv
    params(end+1) = MPCParams(1,'3a(xiv) MassPointTrans');
    params(end).rho_logyP = sqrt(0.9516);
    params(end).sd_logyP = sqrt(0.0434);
    params(end).sd_logyT = sqrt(0.6431);
    params(end).lambdaT = 0.0760;
    
    %----------------------------------------------------------------------
    % PART 3b, QUARTERLY MODEL
    %----------------------------------------------------------------------
    
    % i
    params(end+1) = MPCParams(4,'3b(i) KMPTransf');
    params(end).rho_logyP = 0.9879;
    params(end).sd_logyP = sqrt(0.0109);
    params(end).sd_logyT = sqrt(0.0494);
    
    % iv
    params(end+1) = MPCParams(4,'3b(iv) PersEveryPeriod');
    params(end).rho_logyP = 0.9884;
    params(end).sd_logyP = sqrt(0.0105);
    params(end).sd_logyT = sqrt(1.5298);
    params(end).lambdaT = 0.0813;
    
    %----------------------------------------------------------------------
    % PART 4, Exotic Preferences
    %----------------------------------------------------------------------
    for ifreq = [1 4]
        % temptation
        for itempt = [0.005 0.01 0.05]
            params(end+1) = MPCParams(ifreq,['4 Temptation' num2str(itempt)]);
            params(end).temptation = itempt;
            if itempt == 0.05 && ifreq == 4
                params(end).set_betaH_distance(-1e-5);
            end    
        end
        
        % epstein-zin: vary risk_aver
        for ra = [0.5 0.75 1.5 2 4 8]
            params(end+1) = MPCParams(ifreq,['EZ ra' num2str(ra) ' invies1']);
            params(end).risk_aver = ra;
            params(end).invies = 1;
            params(end).EpsteinZin = 1;
        end
        
        % epstein-zin: vary invies
        for ies = [1/4 1/2 3/4 1.5 2 5]
            params(end+1) = MPCParams(ifreq,['EZ ra1 ies' num2str(ies)]);
            params(end).risk_aver = 1;
            params(end).invies = 1/ies;
            params(end).EpsteinZin = 1;
        end
    end
    
    %----------------------------------------------------------------------
    % ADJUST PARAMETERS FOR FREQUENCY
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);
    
    %----------------------------------------------------------------------
    % SET BETA UPPER BOUND FOR EZ CASES
    %----------------------------------------------------------------------
    % varying risk_aver
    EZ = find([params.EpsteinZin]==1 & [params.invies]==1 & [params.freq]==1);
    for iz = EZ
        params(iz).set_betaH_distance(-3e-2);
    end
    EZ = find([params.EpsteinZin]==1 & [params.invies]==1 & [params.freq]==4);
    for iz = EZ
        params(iz).set_betaH_distance(-6e-3);
    end
    
    
    % --------- annual, varying invies -------------
    freq = 1;
    
    change_betaH = ['EZ ra1 ies' num2str(1/4)];
    params.set_betaH_distance(-3e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(1/2)];
    params.set_betaH_distance(-3e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(3/4)];
    params.set_betaH_distance(-3e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(1.5)];
    params.set_betaH_distance(-2.1e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(2)];
    params.set_betaH_distance(-2.1e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(5)];
    params.set_betaH_distance(-2.1e-2,change_betaH,freq);
    
    % --------- quarterly, varying invies -------------
    freq = 4;
    
    change_betaH = ['EZ ra1 ies' num2str(1/4)];
    params.set_betaH_distance(-1.5e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(1/2)];
    params.set_betaH_distance(-1.5e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(3/4)];
    params.set_betaH_distance(-1e-2,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(1.5)];
    params.set_betaH_distance(-5.5e-3,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(2)];
    params.set_betaH_distance(-5.5e-3,change_betaH,freq);
    
    change_betaH = ['EZ ra1 ies' num2str(5)];
    params.set_betaH_distance(-5.4e-3,change_betaH,freq);
    
    
    %----------------------------------------------------------------------
    % SET BETA UPPER BOUND FOR BETA HETEROGENEITY CASES
    %----------------------------------------------------------------------
    
    % fixed beta heterogeneity cases converge ok
    
    % --------- annual, random ----------
    freq = 1;
    
    % no death
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 NoDeath';
    params.set_betaH_distance(-3e-2,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 NoDeath';
    params.set_betaH_distance(-3e-2,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 NoDeath';
    params.set_betaH_distance(-2e-2,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 NoDeath';
    params.set_betaH_distance(-2e-2,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 NoDeath';
    params.set_betaH_distance(-1.5e-2,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 NoDeath';
    params.set_betaH_distance(-1.5e-2,change_betaH,freq);
        
    % death
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 Death';
    params.set_betaH_distance(-1.3e-2,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 Death';
    params.set_betaH_distance(-9e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 Death';
    params.set_betaH_distance(-7e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 Death';
    params.set_betaH_distance(-4e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 Death';
    params.set_betaH_distance(-1e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 Death';
    params.set_betaH_distance(3e-3,change_betaH,freq);

    % --------- quarterly, random -----------
    freq = 4;
    
    % no death
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 NoDeath';
    params.set_betaH_distance(-1e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 NoDeath';
    params.set_betaH_distance(-1e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 NoDeath';
    params.set_betaH_distance(3e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 NoDeath';
    params.set_betaH_distance(5e-3,change_betaH,freq);
                    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 NoDeath';
    params.set_betaH_distance(4e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 NoDeath';
    params.set_betaH_distance(1.5e-2,change_betaH,freq);
    
    % death
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.02 Death';
    params.set_betaH_distance(-1e-4,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.001 SwitchProb0.1 Death';
    params.set_betaH_distance(-1e-4,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.02 Death';
    params.set_betaH_distance(2e-3,change_betaH,freq);

    change_betaH = '2 RandomBetaHet5 Width0.005 SwitchProb0.1 Death';
    params.set_betaH_distance(5.1e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.02 Death';
    params.set_betaH_distance(4.7e-3,change_betaH,freq);
    
    change_betaH = '2 RandomBetaHet5 Width0.01 SwitchProb0.1 Death';
    params.set_betaH_distance(1.2e-2,change_betaH,freq);

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    if runopts.Server==1 || runopts.Display==0
        [params.Display] = deal(0);
    else
        [params.Display] = deal(1);
    end
    
    if runopts.fast == 1
        params.set_fast();
    end
    
    params = MPCParams.select_by_names(params,selection.names_to_run);
    
    if numel(selection.frequencies) == 1
        params = MPCParams.select_by_freq(params,selection.frequencies);
    end
    
    params.set_index();
    
%     warning('BETA ITERATION IS SET OFF')
%     [params.IterateBeta] = deal(0);
end