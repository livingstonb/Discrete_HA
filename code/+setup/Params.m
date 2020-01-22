classdef Params < handle
    % Usage: params = Params(frequency,name,IncomeProcess)
    % After instantiating a Params object, modify any
    % of the parameters using dot notation. The method
    % 'adjust_if_quarterly' must be called at the end since 
    % it is assumed that the discount factor, returns, etc... 
    % are all entered in annual terms.

    % IncomeProcess is the directory of the income grids,
    % or an empty string to generate the income process
    % in the code.

    properties
        % identifiers
        name;
        index;
        
        % data frequency
        freq;

        % source for income process (file, or empty string for gen in code)
        IncomeProcess ='';
        
        path;

        % computation
        max_iter    = 1e5; % EGP
        tol_iter    = 1.0e-6; % EGP
        Nsim        = 2e5; % For optional simulation
        
        % beta iteration
        maxiterAY   = 50;
        tolAY       = 1e-7;

        % mpc options
        Nmpcsim = 2e5; % Number of draws to compute MPCs
        mpcfrac = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
        
        % wealth statistics options
        epsilon = [0, 0.005, 0.01, 0.02, 0.05, 0.1 0.15]; % fraction of mean ann labor income
        percentiles = [10, 25, 50, 75, 90, 95, 99, 99.9]; % in percent
        
        % decomposition
        abars = [0, 0.01, 0.05];

        % cash on hand / savings grid
        nx          = 500;
        nx_DST      = 400;
        xmax        = 100;
        xgrid_par   = 0.2; %1 for linear, 0 for L-shaped
        borrow_lim  = 0; % negative does not work
        gridspace_min = 0.001; % minimum grid space (0 for no minimum)
        
        % OPTIONS
        MakePlots = 0;
        Simulate = 0;
        MPCs;
        mpcshocks_after_period1;

        % returns
        r = 0.02; % default annual, adjusted if frequency = 4;
        R;
        annuities = false;
        
        % demographics
        dieprob = 1/50; % default annual, adjusted if frequency = 4;
        
        % preferences
        EpsteinZin  = 0;
        invies      = 2.5; % only relevant for Epstein-Zin
        risk_aver   = 1;
    	beta0       = 0.98;
    	temptation  = 0;
    	betaL       = 0.80;
        betaH0      = - 1e-3; % adjustment factor if not using theoretical upper bound
        betaH; % theoretical upper bound if not adjusting
        
        % warm glow bequests: bequest weight = 0 is accidental
        bequest_weight  = 0;
        bequest_curv    = 1;
        bequest_luxury  = 0.01;
        Bequests = 1; % 1 for wealth left as bequest, 0 for disappears
        
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
        ResetIncomeUponDeath = 1;
        
        % government
        labtaxlow       = 0; %proportional tax
        labtaxhigh      = 0; %additional tax on incomes above threshold
        labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in
        savtax          = 0; %0.0001;  %tax rate on savings
        savtaxthresh    = 0; %multiple of mean gross labor income
        lumptransfer   = 0;

        % discount factor shocks
        nbeta = 1;
        beta_dist = 1; % either 1 for equal prob in all, or vector summing to 1
        betawidth = 0.005;
        betaswitch = 0;
        beta_grid_forced = []; % overrides all other beta values if used

        % used for different het cases, need to recode this
        nb = 1;

        % IES shocks
        IESswitch = 0;
        
        % computation
    	Tsim        = 400; % Simulation

        % beta iteration
        IterateBeta;
    	targetAY    = 3.5; 

    	% mpc shock sizes
    	shocks = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
    end

    methods
        function obj = Params(frequency,name,IncomeProcess)
        	% create params object
            obj.name = name;
            obj.freq = frequency;
            obj.IncomeProcess = IncomeProcess;
            obj.R = 1 + obj.r;
            
            if frequency == 1
                obj.sd_logyT = sqrt(0.0494);
                obj.sd_logyP = sqrt(0.0422);
                obj.rho_logyP =0.9525;
            elseif frequency == 4
                % use quarterly_a if quarterly & no IncomeProcess is given
                obj.sd_logyT = sqrt(0.2087);
                obj.sd_logyP = sqrt(0.01080);
                obj.rho_logyP = 0.9881;
            else
                error('Frequency must be 1 or 4')
            end             
        end
        
        function obj = annuities_on(obj)
            % Wait until after frequency adjustments to change other
            % variables
            if numel(obj) > 1
                error('Turn on annuities for single specification at a time only')
            end
            obj.annuities = 1;
        end

        function obj = set_run_parameters(obj,runopts)
        	% use fields in runopts to set values in Params object

            % fast option
            if runopts.fast == 1
                [obj.nx_DST] = deal(10);
                [obj.nx] = deal(10);
                [obj.Nmpcsim] = deal(1e2);
                [obj.nyT] = deal(3);
                [obj.nyP] = deal(3);
                [obj.Tsim] = deal(100);
            end
            
            % iterate option
            [obj.IterateBeta] = deal(runopts.IterateBeta);

            % simulate option
            [obj.Simulate] = deal(runopts.Simulate);

            % compute mpcs
            [obj.MPCs] = deal(runopts.MPCs);

            % compute mpcs for is > 1?
            [obj.mpcshocks_after_period1] = deal(runopts.mpcshocks_after_period1);
            
            [obj.path] = deal(runopts.path);
        end
        
        function obj = set_index(obj)
        	% reset index to count from 1 to numel(params)
            ind = num2cell(1:numel(obj));
            [obj.index] = deal(ind{:});
        end

        function obj = set_grid(obj,nx,nxlong,curv)
        	% convenient way to set grid for grid tests
            obj.nx = nx;
            obj.nxlong = nxlong;
            obj.xgrid_par = curv;
        end
    end
    
    methods (Static)
        
        function objs = adjust_if_quarterly(objs)
        	% adjusts relevant parameters such as r to a quarterly frequency,
        	% i.e. r = r / 4
            % must be called after setting all parameterizations

            for io = 1:numel(objs)
                objs(io).R = (1+objs(io).r).^(1/objs(io).freq);
                objs(io).r = objs(io).R - 1;

                objs(io).savtax = objs(io).savtax/objs(io).freq;
                objs(io).Tsim = objs(io).Tsim * objs(io).freq; % Increase simulation time if quarterly
                objs(io).beta0 = objs(io).beta0^(1/objs(io).freq);
                objs(io).dieprob = 1 - (1-objs(io).dieprob)^(1/objs(io).freq);
                objs(io).betaswitch = 1 - (1-objs(io).betaswitch)^(1/objs(io).freq);

                objs(io).betaL = objs(io).betaL^(1/objs(io).freq);
                
                objs(io).betaH = 1./((max(objs(io).R))*(1-objs(io).dieprob));
                objs(io).betaH = objs(io).betaH + objs(io).betaH0;
	                
                if objs(io).annuities == true
                    objs(io).Bequests = 0;
                    objs(io).r = objs(io).r + objs(io).dieprob;
                    objs(io).R = 1 + objs(io).r;
                end

                objs(io).lumptransfer = objs(io).lumptransfer / objs(io).freq;
            end
        end
        
        function objs = select_by_number(objs,number)
            objs = objs([objs.index]==number);
        end
        
        function objs = select_by_names(objs,names_to_run)
            % discards all experiments with names not included in the
            % cell array 'names_to_run'
            if ~isempty(names_to_run)
                % Indices of selected names within params
                params_to_run = ismember({objs.name},names_to_run);
                if sum(ismember(names_to_run,{objs.name})) < numel(names_to_run)
                    error('Some of the entries in names_to_run are invalid')
                else
                    objs = objs(params_to_run);
                end
            else
                return
            end
        end
    end

end
