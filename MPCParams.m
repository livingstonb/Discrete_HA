classdef MPCParams < handle
    % usage: params = MPCParams(frequency)
    
    properties % default properties
        
        % common properties
        filesuffix;
        Fast;
        Server;
        
        % identifiers
        name;
        index;
        
        % data frequency
        freq;
        
        % returns
        r = 0.02; % default annual, adjusted if frequency = 4;
        R;
        
        % demographics
        dieprob = 1/50; % default annual, adjusted if frequency = 4;
        
        % preferences
        EpsteinZin  = 0;
        invies      = 2.5;
        risk_aver   = 1;
    	beta0       = 0.98;
    	temptation  = 0;
    	betaL       = 0.80;
        betaH;
        
        % warm glow bequests: bequest weight = 0 is accidental
        bequest_weight  = 0;
        bequest_curv    = 1;
        bequest_luxury  = 0.01;
        Bequests = 1; % 1 for wealth left as bequest, 0 for disappears
        Annuities = 0; % Automatically turns off bequests if set to 1
        
        % source for income process (file, or empty string for gen in code)
        IncomeProcess ='';
        
        % income risk: AR(1) + IID in logs
    	nyT	= 11;
        yTContinuous = 0; % only works for simulation
        sd_logyT;
        lambdaT = 1;
        nyP = 11; % persistent component
        sd_logyP;
        rho_logyP;
        nyF = 1;
        sd_logyF = 0;
        
        % cash on hand / savings grid
    	nx          = 100;
    	xmax        = 1000;
    	xgrid_par   = 1/3; %1 for linear, 0 for L-shaped
    	borrow_lim  = 0; % negative does not work
        
        %government
        labtaxlow       = 0; %proportional tax
        labtaxhigh      = 0; %additional tax on incomes above threshold
        labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in
        savtax          = 0; %0.0001;  %tax rate on savings
        savtaxthresh    = 0; %multiple of mean gross labor income

        %discount factor shocks
        nb = 1;
        betawidth = 0.005;
        betaswitch = 0;
        
        % computation
    	max_iter    = 1e5; % EGP
    	tol_iter    = 1.0e-6; % EGP
    	Nsim        = 100000; % For optional simulation
    	Tsim        = 200; % Simulaation
    	nxlong      = 750;

        % beta iteration
    	targetAY    = 3.5;
    	maxiterAY   = 50;
    	tolAY       = 1e-5;
        
        % mpc options
        Nmpcsim = 1e6; % Number of draws to compute MPCs
        mpcfrac = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
        
        % wealth statistics options
    	epsilon = [0, 0.005, 0.01, 0.02, 0.05, 0.1]; % fraction of mean ann labor income
    	percentiles = [10, 25, 50, 75, 90, 95, 99, 99.9]; % in percent
        
        % decomposition
        abars = [0, 0.01, 0.05];
        
        % OPTIONS
    	IterateBeta = 1;
        Display = 0;
        MakePlots = 0;
        Simulate = 0;
        
    end

    methods
        function obj = MPCParams(frequency,name)
            obj.name = name;
            
            % Adjust for frequency
            if frequency == 1
                obj.freq = 1;
                obj.R = 1 + obj.r;
                
                obj.sd_logyT = 0.0497;
                obj.sd_logyP = 0.0422;
                obj.rho_logyP =0.9525;
            elseif frequency == 4;
                obj.freq = 4;
            
                obj.R = (1+obj.r)^(1/obj.freq);
                obj.r = obj.R - 1;

                obj.savtax        = obj.savtax/obj.freq;
                obj.Tsim          = obj.Tsim * obj.freq; % Increase simulation time if quarterly
                obj.beta0         = obj.beta0^(1/obj.freq);
                obj.dieprob       = 1 - (1-obj.dieprob)^(1/obj.freq);
                obj.betaswitch    = 1 - (1-obj.betaswitch)^(1/obj.freq);
                obj.betaL         = obj.betaL^(1/obj.freq);
                obj.betaH         = 1/((obj.R)*(1-obj.dieprob));
                
                obj.sd_logyT = 0.02087;
                obj.sd_logyP = 0.0108;
                obj.rho_logyP =0.9881;
            else
                error('Frequency must be 1 or 4')
            end
        end
        
        function obj = set_fast(obj)
            [obj.nxlong] = deal(20);
            [obj.nx] = deal(15);
            [obj.Nmpcsim] = deal(1e2);
            [obj.nyT] = deal(3);
            [obj.nyP] = deal(3);
        end
        
        function obj = set_index(obj)
            % Should be called only after selecting which runs to drop
            ind = num2cell(1:numel(obj));
            [obj.index] = deal(ind{:});
        end
        
        function obj = set_filesuffix(obj,filesuffix)
            % filesuffix should be of the form '_quarterlyA'
            [obj.filesuffix] = deal(filesuffix);
        end
      
    end
    
    methods (Static)  
        
        function objs = select_by_names(objs,names_to_run)
            % Choose parameterizations based on name
            if isequal(class(names_to_run),'cell')
                if isempty(names_to_run)
                    return 
                else
                    % Indices of selected names within params
                    params_to_run = ismember({objs.name},names_to_run);
                    objs = objs(params_to_run);
                end
            else
                error('names_to_run must be a cell array')
            end
        end
        
        function objs = select_by_freq(objs,freq)
            % Choose parameterizations based on frequency
            if freq==1 || freq==4
                objs = objs([objs.freq]==freq);
            else
                error('Must select freq = 1 or freq = 4')
            end
        end 
        
        function objs = add(objs,newparams)
            objs(end+1) = [];
        end
        
    end

end