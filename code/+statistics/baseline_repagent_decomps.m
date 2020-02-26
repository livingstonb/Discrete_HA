function [decomp_wrt_baseline, decomp_wrt_repagent] = baseline_and_repagent_decomps(params, results, return_nans)
    % decomp_wrt_baseline is relative to baseline
    % decomp_wrt_repagent is relative to representative agent model
    
    % computed with respect to a shock of 0.01 * mean ann gross income

    % Construct agrid based off params(1)
    agrid = results(1).direct.agrid;
    agrid_short = agrid;
    agrid = repmat(agrid, params(1).nyP*params(1).nyF*params(1).nb, 1);

    % Initialize with NaNs
    for ip = 1:numel(params)
        decomp_wrt_baseline(ip).mpc1_Em1_less_Em0 = NaN;
        decomp_wrt_baseline(ip).mpc1_term1 = NaN;
        decomp_wrt_baseline(ip).mpc1_term2 = NaN;
        decomp_wrt_baseline(ip).mpc1_term3 = NaN;
        decomp_wrt_repagent(ip).mpc1_Em1_less_mRA = NaN;
        decomp_wrt_repagent(ip).mpc1_term1 = NaN;
        decomp_wrt_repagent(ip).mpc1_term2 = NaN;
        decomp_wrt_repagent(ip).mpc1_term3 = NaN;
        
        decomp_wrt_baseline(ip).mpc4_Em1_less_Em0 = NaN;
        decomp_wrt_baseline(ip).mpc4_term1 = NaN;
        decomp_wrt_baseline(ip).mpc4_term2 = NaN;
        decomp_wrt_baseline(ip).mpc4_term3 = NaN;
        decomp_wrt_repagent(ip).mpc4_Em1_less_mRA = NaN;
        decomp_wrt_repagent(ip).mpc4_term1 = NaN;
        decomp_wrt_repagent(ip).mpc4_term2 = NaN;
        decomp_wrt_repagent(ip).mpc4_term3 = NaN;

        decomp_wrt_baseline(ip).mpc1_term2a = NaN(numel(params(1).abars),1);
        decomp_wrt_baseline(ip).mpc1_term2b = NaN(numel(params(1).abars),1);
        decomp_wrt_baseline(ip).mpc4_term2a = NaN(numel(params(1).abars),1);
        decomp_wrt_baseline(ip).mpc4_term2b = NaN(numel(params(1).abars),1);

        decomp_wrt_baseline(ip).decomp_error = NaN;
        decomp_wrt_repagent(ip).decomp_error = NaN;
        
        if return_nans || (~params(ip).MPCs)
            continue
        end

        % find index of baseline within 'p'
        if params(ip).freq == 1
            baseind = find(ismember({params.name}, {'baseline_A'}));
        else
            baseind = find(ismember({params.name}, {'baseline_Q'}));
        end
        
        skip_baseline_decomp = (numel(params) == 1) || isempty(baseind);

        if ~skip_baseline_decomp
            %% --------------------------------------------------------------
    		% BASELINE DISTRIBUTION AND MPCs OVER ASSET GRID
    		% ---------------------------------------------------------------
    		dim0 = params(baseind).nyF * params(baseind).nyP * params(baseind).nb;

    		% distribution
            g0 = results(baseind).direct.adist(:); 
            
            % one-period mpcs
            m0 = results(baseind).direct.mpcs(5).mpcs_1_t{1};

            % dist and one-period mpcs condensed to asset grid
            [g0_a, m0_a] = get_functions_of_assets(baseind, params, g0, m0, dim0);

            % four-period mpcs
            m0_4 = 0;
            for t = 1:4
                m0_4 = m0_4 + results(baseind).direct.mpcs(5).mpcs_1_t{t};
            end

            % four-period mpcs condensed to asset grid
            [~,m0_4_a] = get_functions_of_assets(baseind, params, g0, m0_4, dim0);
        end

        %% --------------------------------------------------------------
		% ALTERNATIVE MODEL DISTRIBUTION AND MPCs OVER ASSET GRID
		% ---------------------------------------------------------------
        if params(ip).nb > 1
            continue
        end
        
        dim1 = params(ip).nyF * params(ip).nyP * params(ip).nb;

        % distribution
        g1 = results(ip).direct.adist(:);

        % mpcs
        m1 = results(ip).direct.mpcs(5).mpcs_1_t{1};       

        % distribution and mpcs condensed to asset grid
        [g1_a, m1_a] = get_functions_of_assets(ip, params, g1, m1, dim1);

        % 4-period mpcs
        if params(ip).freq == 4
            m1_4 = 0;
            for i = 1:4
                % model MPCs
                m1_it = results(ip).direct.mpcs(5).mpcs_1_t{i};
                m1_4 = m1_4 + m1_it;
            end

            [~, m1_4_a] = get_functions_of_assets(ip, params, g1, m1_4, dim1);
        end

        if ~skip_baseline_decomp
            %% --------------------------------------------------------------
    		% DECOMP WITH RESPECT TO BASELINE
    		% ---------------------------------------------------------------
            decomp_wrt_baseline(ip).mpc1_Em1_less_Em0 = results(ip).direct.mpcs(5).avg_s_t(1,1) ...
                            - results(baseind).direct.mpcs(5).avg_s_t(1,1);
            decomp_wrt_baseline(ip).mpc1_term1 = g0_a' * (m1_a - m0_a); % Effect of MPC function
            decomp_wrt_baseline(ip).mpc1_term2 = m0_a' * (g1_a - g0_a); % Effect of distribution
            decomp_wrt_baseline(ip).mpc1_term3 = (m1_a - m0_a)' * (g1_a - g0_a); % Interaction

            % Decomposition of distribution effect
            % create interpolants from assets to integrals of m0g0, m0g1
            m0g0 = m0_a .* g0_a;
            m0g0_interp = griddedInterpolant(agrid_short,cumsum(m0g0),'linear');
            m0g1 = m0_a .* g1_a;
            m0g1_interp = griddedInterpolant(agrid_short,cumsum(m0g1),'linear');

            for ia = 1:numel(params(ip).abars)
                abar = params(ip).abars(ia); % HtM threshold
                if abar == 0
    	            idx = agrid_short <= abar;
    	            % HtM households
    	            decomp_wrt_baseline(ip).mpc1_term2a(ia) = m0_a(idx)' * (g1_a(idx) - g0_a(idx));
    	            % Non-HtM households
    	            decomp_wrt_baseline(ip).mpc1_term2b(ia) = m0_a(~idx)' * (g1_a(~idx) - g0_a(~idx));
    	        else
    	        	decomp_wrt_baseline(ip).mpc1_term2a(ia) = m0g1_interp(abar) - m0g0_interp(abar);
    	        	decomp_wrt_baseline(ip).mpc1_term2b(ia) = (m0_a'*g1_a - m0g1_interp(abar)) ...
    	        									- (m0_a'*g0_a - m0g0_interp(abar));
    	        end
            end

            %% 4-Period MPC decomp_wrt_baseline (decomp with respect to baseline)
            % Only used for quarterly models
            if params(ip).freq == 4
                % create interpolants from assets to m0g0, m0g1
                m0g0_4 = m0_4_a .* g0_a;
                m0g0_interp4 = griddedInterpolant(agrid_short,cumsum(m0g0_4),'linear');
                m0g1_4 = m0_4_a .* g1_a;
                m0g1_interp4 = griddedInterpolant(agrid_short,cumsum(m0g1_4),'linear');

                decomp_wrt_baseline(ip).mpc4_Em1_less_Em0 = results(ip).direct.mpcs(5).avg_1_1to4 ...
                            - results(baseind).direct.mpcs(5).avg_1_1to4;
                decomp_wrt_baseline(ip).mpc4_term1 = g0_a' * (m1_4_a - m0_4_a);
                decomp_wrt_baseline(ip).mpc4_term2 = m0_4_a' * (g1_a - g0_a);
                decomp_wrt_baseline(ip).mpc4_term3 = (m1_4_a - m0_4_a)' * (g1_a - g0_a);

                for ia = 1:numel(params(ip).abars)
                    abar = params(ip).abars(ia);
                    if abar == 0
                        idx = agrid_short <= abar;
                        decomp_wrt_baseline(ip).mpc4_term2a(ia) = m0_4_a(idx)' * (g1_a(idx) - g0_a(idx));
                        decomp_wrt_baseline(ip).mpc4_term2b(ia) = m0_4_a(~idx)' * (g1_a(~idx) - g0_a(~idx));
                    else
                        decomp_wrt_baseline(ip).mpc4_term2a(ia) = m0g1_interp4(abar) - m0g0_interp4(abar);
                        decomp_wrt_baseline(ip).mpc4_term2b(ia) = (m0_4_a'*g1_a - m0g1_interp4(abar)) ...
                                                    - (m0_4_a'*g0_a - m0g0_interp4(abar));
                    end
                end
            end
        end

        %% --------------------------------------------------------------
		% DECOMP WITH RESPECT TO REP AGENT MODEL
		% ---------------------------------------------------------------
        if params(ip).nb == 1
            % RA MPC
            temp = (1-params(ip).dieprob) * results(ip).direct.beta * params(ip).R;
            m0_RA = params(ip).R * temp ^ (-1/params(ip).risk_aver) - 1;

            % To get E[MPC|a=3.5], interpolate
            m1_a_interp = griddedInterpolant(agrid_short,m1_a,'linear');
            m1_atmean = m1_a_interp(3.5);

            % Mean model MPC - RA MPC
            decomp_wrt_repagent(ip).mpc1_Em1_less_mRA = results(ip).direct.mpcs(5).avg_s_t(1,1) - m0_RA;

            decomp_wrt_repagent(ip).mpc1_term1 = m1_atmean - m0_RA; % Effect of MPC function
            decomp_wrt_repagent(ip).mpc1_term2 = 0; % Effect of distribution
            decomp_wrt_repagent(ip).mpc1_term3 = (m1_a - m0_RA)' * g1_a - (m1_atmean - m0_RA); % Interaction
        end
    end
end

function [pmf_assets, mpcs_assets] = get_functions_of_assets(baseind, params, pmf, mpcs, dim)
	% this function outputs pmf(a) and E[mpc|a], functions of assets,
	% given functions over the entire state space
	
	pmf_wide = reshape(pmf, [params(baseind).nx_DST dim]);
	pmf_assets = sum(pmf_wide,2);
	pmf_small = pmf_assets < 1e-9;

	mpcs_wide = reshape(mpcs,[params(baseind).nx_DST dim]);
	mpcs_assets = mpcs_wide .* pmf_wide ./ pmf_assets;
	mpcs_assets = sum(mpcs_assets,2);

	% at asset values where P(a) is too small to divide accurately,
	% replace expected value with arithmetic mean of mpcs
	mpcs_assets(pmf_small) = mean(mpcs_wide(pmf_small,:),2);
end