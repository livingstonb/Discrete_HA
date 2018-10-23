classdef MPCParams < handle
    % usage: params = MPCParams(frequency)
    
    properties
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
        betaH0;
        betaH;
        
        % warm glow bequests: bequest weight = 0 is accidental
        bequest_weight  = 0;
        bequest_curv    = 1;
        bequest_luxury  = 0.01;
        Bequests = 1; % 1 for wealth left as bequest, 0 for disappears
        
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
        Display;
        MakePlots = 0;
        Simulate = 0;
        
    end

    methods
        function obj = MPCParams(frequency,name)
            obj.name = name;
            obj.freq = frequency;
            obj.R = 1 + obj.r;
            
            if frequency == 1
                obj.sd_logyT = sqrt(0.0497);
                obj.sd_logyP = sqrt(0.0422);
                obj.rho_logyP =0.9525;
            elseif frequency == 4
                obj.sd_logyT = sqrt(0.02087);
                obj.sd_logyP = sqrt(0.0108);
                obj.rho_logyP = 0.9881;
            else
                error('Frequency must be 1 or 4')
            end             
        end
        
        function obj = set_betaH_distance(obj,val,name,freq)
            if nargin == 4
                change_ind = find(ismember({obj.name},name) & [obj.freq]==freq);
                obj(change_ind).betaH = obj(change_ind).betaH0 + val;
            elseif nargin == 2
                if numel(obj) == 1
                    obj.betaH = obj.betaH0 + val;
                else
                    error('cannot pass array to method unless passing name,freq')
                end
            else 
                error('must pass 1 or 3 arguments')
            end
        end
        
        function obj = annuities_on(obj)
            % Turn off bequests
            if numel(obj) > 1
                error('This method only adjusts one parameterization at a time')
            else
                obj.Bequests = 0;
                obj.r = obj.r + obj.dieprob;
                obj.R = 1 + obj.r;
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
      
    end
    
    methods (Static)
        
        function objs = adjust_if_quarterly(objs)
            for io = 1:numel(objs)
                objs(io).R = (1+objs(io).r)^(1/objs(io).freq);
                objs(io).r = objs(io).R - 1;

                objs(io).savtax = objs(io).savtax/objs(io).freq;
                objs(io).Tsim = objs(io).Tsim * objs(io).freq; % Increase simulation time if quarterly
                objs(io).beta0 = objs(io).beta0^(1/objs(io).freq);
                objs(io).dieprob = 1 - (1-objs(io).dieprob)^(1/objs(io).freq);
                objs(io).betaswitch = 1 - (1-objs(io).betaswitch)^(1/objs(io).freq);
                objs(io).betaL = objs(io).betaL^(1/objs(io).freq);
                
                objs(io).betaH0 = 1/((objs(io).R)*(1-objs(io).dieprob));
                objs(io).betaH = objs(io).betaH0 - 1e-3;
            end
        end
        
        function objs = select_by_names(objs,names_to_run)
            % Choose parameterizations based on name
            if ~isempty(names_to_run)
                % Indices of selected names within params
                params_to_run = ismember({objs.name},names_to_run);
                if sum(ismember(names_to_run,{objs.name})) < numel(names_to_run)
                    disp('Valid specification names include:')
                    % valid_names = {objs.name};
                    fprintf('%s \n',valid_names{:})
                    error('Some of the entries in names_to_run are invalid')
                else
                    objs = objs(params_to_run);
                end
            else
                return
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
        
    end

end