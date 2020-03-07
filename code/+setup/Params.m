classdef Params < handle
    % Usage: params = Params(frequency,name,IncomeProcess)
    % After instantiating a Params object, modify any
    % of the parameters using dot notation. The method
    % 'adjust_if_quarterly' must be called at the end since 
    % it is assumed that the discount factor, returns, etc... 
    % are all entered in annual terms.
    %
    % IncomeProcess is the directory of the income grids,
    % or an empty string to generate the income process
    % in the code.
    %
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    properties
        % identifiers
        name;
        index;
        
        % data frequency
        freq;

        % source for income process (file, or empty string for gen in code)
        IncomeProcess = '';
        
        path;

        % mean annual income dollar interpretation
        annual_inc_dollars = 72000;
        convert_to_dollars;
        convert_from_dollars;

        % computation
        max_iter = 1e5; % EGP
        tol_iter = 1.0e-6; % EGP
        Nsim = 1e5; % For optional simulation

        % mpc options
        shocks = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
        shocks_labels; % labels
        Nmpcsim = 2e5; % Number of draws to compute MPCs
        
        % wealth statistics options
        epsilon = [0, 0.005, 0.01, 0.02, 0.05, 0.1 0.15]; % fraction of mean ann labor income
        percentiles = [10, 25, 50, 75, 90, 95, 99, 99.9]; % in percent
        
        % decomposition
        abars = [0, 0.01, 0.05];

        % cash on hand / savings grid
        nx = 250;
        nx_neg = 0;
        nx_DST = 250;
        nx_neg_DST = 0;
        xmax = 50;
        xgrid_par = 0.2; %1 for linear, 0 for L-shaped
        xgrid_par_neg = 0.4;
        borrow_lim = 0;
        nbl_adjustment = 0.99;
        gridspace_min = 0; % minimum grid space (0 for no minimum)
        alternate_gcurv = false;
        
        % OPTIONS
        MakePlots = 0;
        Simulate = 0;
        MPCs = 0;
        MPCs_news = 0;
        MPCs_loan_and_loss = 0;
        DeterministicMPCs = 0;
        outdir = '';
        savematpath = '';
        SaveOutput = false;

        % returns
        r = 0.02; % default annual, adjusted if frequency = 4;
        R;
        annuities = false;
        
        % demographics
        dieprob = 1/50; % default annual, adjusted if frequency = 4;
        
        % preferences
        EpsteinZin = 0;
        invies = 2.5; % only relevant for Epstein-Zin
        risk_aver = 1;
    	beta0 = 0.98;
    	temptation = 0;
    	betaL = 0.80;
        betaH0 = -1e-3; % adjustment factor if not using theoretical upper bound
        betaH; % theoretical upper bound if not adjusting
        
        % warm glow bequests: bequest weight = 0 is accidental
        bequest_weight = 0;
        bequest_curv = 1;
        bequest_luxury = 0.01;
        Bequests = true; % 1 for wealth left as bequest, 0 for disappears
        
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
        ResetIncomeUponDeath = true;
        
        % government
        labtaxlow = 0; %proportional tax
        labtaxhigh = 0; %additional tax on incomes above threshold
        labtaxthreshpc = 0.99; %percentile of earnings distribution where high tax rate kicks in
        savtax = 0; %0.0001;  %tax rate on savings
        savtaxthresh = 0; %multiple of mean gross labor income
        compute_savtax;
        lumptransfer = 0;

        % discount factor shocks
        nbeta = 1;
        beta_dist = 1; % either 1 for equal prob in all, or vector summing to 1
        betawidth = 0.005;
        beta_grid_forced = []; % overrides all other beta values if used

        % used for different het cases
        nb = 1;

        % probability of switch in z-dimension
        prob_zswitch = 0;
        zdist_forced;
        
        % computation
    	Tsim = 400; % Simulation

        % calibration
        calibrate;
        calibrate_maxiter = 60;
        calibrate_tol = 1e-6;
        calibrator;
        x0_calibration;
        target_value = 3.5;

        % other, unspecified option
        other;
    end

    methods
        function obj = Params(frequency, name, IncomeProcess)
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

        function obj = set_run_parameters(obj, runopts)
        	% use fields in runopts to set values in Params object

            % fast option
            if runopts.fast
                obj.nx_DST = 10;
                obj.nx = 11;
                obj.Nmpcsim = 1e2;
                obj.nyT = 5;
                obj.nyP = 3;
                obj.Tsim = 100;
            end

            obj.MakePlots = runopts.MakePlots;

            obj.outdir = runopts.outdir;

            obj.savematpath = runopts.savematpath;
            
            % iterate option
            obj.calibrate = runopts.calibrate;

            % simulate option
            obj.Simulate = runopts.Simulate;

            % compute mpcs
            obj.MPCs = runopts.MPCs;

            % compute mpcs for is > 1?
            obj.MPCs_news = runopts.MPCs_news;

            obj.MPCs_loan_and_loss = runopts.MPCs_loan_and_loss;
            obj.DeterministicMPCs = runopts.DeterministicMPCs;
            obj.path = runopts.path;
        end
        
        function obj = set_index(obj)
        	% reset index to count from 1 to numel(params)
            ind = num2cell(1:numel(obj));
            [obj.index] = deal(ind{:});
        end

        function obj = set(obj, field, new_val, quiet)
            % Sets the value of a parameter.
            %
            % Inputs
            % ------
            %
            % field : A string containing the parameter
            %   name.
            %
            % new_val : The desired value of the parameter.
            %
            % quiet : An optional argument that, when it
            %   evaluates to true, suppresses printing
            %   to the screen.

            field = char(field);
            if ~isprop(obj, field)
                error("Requested field is not an attribute of Params.");
            end

            obj.(field) = new_val;

            if ~exist('quiet', 'var')
                quiet = false;
            end

            if ~quiet
                disp(strcat(field, sprintf(" has been reset to %.9f", new_val)));
            end
        end

        function [matches, indices] = return_indices(obj, names)
            matches = ismember(names, {obj.name});
            indices = find(ismember({obj.name}, names));
        end

        function make_adjustments(obj)
            obj.make_frequency_adjustments();
            obj.make_other_adjustments();
        end

        function make_frequency_adjustments(obj)
            % adjusts relevant parameters such as r to a quarterly frequency,
            % i.e. r = r / 4
            % must be called after setting all parameterizations

            obj.R = (1 + obj.r) .^ (1 / obj.freq);
            obj.r = obj.R - 1;

            obj.savtax = obj.savtax / obj.freq;
            obj.Tsim = obj.Tsim * obj.freq; % Increase simulation time if quarterly
            obj.beta0 = obj.beta0^(1 / obj.freq);
            obj.dieprob = 1 - (1 - obj.dieprob) ^ (1 / obj.freq);
            obj.prob_zswitch = 1 - (1 - obj.prob_zswitch) ^ (1 / obj.freq);

            obj.betaL = obj.betaL^(1 / obj.freq);

            obj.betaH = 1 ./ (max(obj.R) * (1-obj.dieprob));
            obj.betaH = obj.betaH + obj.betaH0;

            obj.lumptransfer = obj.lumptransfer / obj.freq;
        end

        function make_other_adjustments(obj)
            if obj.annuities
                obj.Bequests = false;
                obj.r = obj.r + obj.dieprob;
                obj.R = 1 + obj.r;
            end

            obj.nbeta = max(obj.nbeta, numel(obj.beta_grid_forced));
            obj.compute_savtax =...
                @(sav) obj.savtax * max(sav - obj.savtaxthresh, 0);
            obj.convert_to_dollars = @(num) num * obj.annual_inc_dollars;
            obj.convert_from_dollars = @(dollars) dollars / obj.annual_inc_dollars;

            if obj.EpsteinZin
                obj.DeterministicMPCs = false;
            end

            if isempty(obj.shocks_labels)
                obj.shocks_labels = {};
                for ishock = 1:6
                    obj.shocks_labels{ishock} = sprintf('%g', obj.shocks(ishock));
                end
            end
        end
    end
    
    methods (Static)
        function objs = select_by_names(objs, names_to_run)
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
