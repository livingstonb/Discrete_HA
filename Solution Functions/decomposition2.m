function [decomp2,decomp3] = decomposition2(params,results)
    % decomp2 is relative to baseline
    % decomp3 is relative to RA model
    
    % computed with respect to a shock of 0.01 * mean ann gross income

    % Construct agrid based off p(1)
    agrid = linspace(0,1,params(1).nx_KFE)';
    agrid = agrid.^(1/params(1).xgrid_par);
    agrid = params(1).borrow_lim + (params(1).xmax - params(1).borrow_lim) * agrid;

    % Force grid spacing >= gridspace_min near 0
    for ia = 1:params(1).nx_KFE-1
        if agrid(ia+1) - agrid(ia) < params(1).gridspace_min
            agrid(ia+1) = agrid(ia) + params(1).gridspace_min;
        else
            break
        end
    end

    agrid_short = agrid;
    agrid = repmat(agrid,params(1).nyP*params(1).nyF*params(1).nb,1);

    % Initialize with NaNs
    for ip = 1:numel(params)
        decomp2(ip).mpc1_Em1_less_Em0 = NaN;
        decomp2(ip).mpc1_term1 = NaN;
        decomp2(ip).mpc1_term2 = NaN;
        decomp2(ip).mpc1_term3 = NaN;
        decomp3(ip).mpc1_Em1_less_mRA = NaN;
        decomp3(ip).mpc1_term1 = NaN;
        decomp3(ip).mpc1_term2 = NaN;
        decomp3(ip).mpc1_term3 = NaN;
        
        decomp2(ip).mpc4_Em1_less_Em0 = NaN;
        decomp2(ip).mpc4_term1 = NaN;
        decomp2(ip).mpc4_term2 = NaN;
        decomp2(ip).mpc4_term3 = NaN;
        decomp3(ip).mpc4_Em1_less_mRA = NaN;
        decomp3(ip).mpc4_term1 = NaN;
        decomp3(ip).mpc4_term2 = NaN;
        decomp3(ip).mpc4_term3 = NaN;

        decomp2(ip).mpc1_term2a = NaN(numel(params(1).abars),1);
        decomp2(ip).mpc1_term2b = NaN(numel(params(1).abars),1);
        decomp2(ip).mpc4_term2a = NaN(numel(params(1).abars),1);
        decomp2(ip).mpc4_term2b = NaN(numel(params(1).abars),1);

        decomp2(ip).decomp_error = NaN;
        decomp3(ip).decomp_error = NaN;

        % find index of baseline within 'p'
        if params(ip).freq == 1
            baseind = find(ismember({params.name},{'baseline_A'}));
        else
            baseind = find(ismember({params.name},{'baseline_Q'}));
        end
        
        %% baseline model
        g0 = results(baseind).direct.adist(:); 

        % g0_wide rows indexed by assets, columns by yF, yP, beta combination
        dim0 = params(baseind).nyF * params(1).nyP * params(baseind).nb;
        g0_wide = reshape(g0,[params(baseind).nx_KFE dim0]);

        % P(a), vector
        g0_a = sum(g0_wide,2);

        % Avoid dividing by very small numbers
        g0small = g0_a < 1e-9;
        g0_a_nonzero = g0_a .* (~g0small) + 2 .* g0small;

        % baseline MPCs
        m0 = results(baseind).direct.mpcs(5).mpcs_1_t{1};       
        m0_wide = reshape(m0,[params(baseind).nx_KFE dim0]);
        m0_a = m0_wide .* g0_wide ./ g0_a_nonzero;
        m0_a = sum(m0_a,2);

        % P(a) too small to divide by accurately, use arithmetic mean
        m0_a(g0small) = mean(m0_wide(g0small,:),2);

        m0_4 = 0;
        for t = 1:4
            m0_4 = m0_4 + results(baseind).direct.mpcs(5).mpcs_1_t{t};
        end
        m0_4_wide = reshape(m0_4,[params(baseind).nx_KFE dim0]);
        m0_4_a = m0_4_wide .* g0_wide ./ g0_a_nonzero;
        m0_4_a = sum(m0_4_a,2);

        % P(a) too small to divide by accurately, use arithmetic mean
        m0_4_a(g0small) = mean(m0_4_wide(g0small,:),2);
        
        %% Model 'ip'
        g1 = results(ip).direct.adist(:);
        
        dim1 = params(ip).nyF * params(ip).nyP * params(ip).nb;
        
        % rows indexed by assets, columns by yF, yP, beta combination
        g1_wide = reshape(g1,[params(ip).nx_KFE dim1]);
        
        % P(a), vector
        g1_a = sum(g1_wide,2);
        
        % Avoid dividing by very small numbers
        g1small = g1_a < 1e-9;
        g1_a_nonzero = g1_a .* (~g1small) + 2 .* g1small;

        
        %% 1-Period MPC decomposition2 (decomp with respect to baseline)
        
        % model MPCs
        m1 = results(ip).direct.mpcs(5).mpcs_1_t{1};       
        m1_wide = reshape(m1,[params(ip).nx_KFE dim1]);
        m1_a = m1_wide .* g1_wide ./ g1_a_nonzero;
        m1_a = sum(m1_a,2);
        
        % P(a) too small to divide by accurately, use arithmetic mean
        m1_a(g1small) = mean(m1_wide(g1small,:),2);

        % Mean model MPC  - mean baseline MPC
        decomp2(ip).mpc1_Em1_less_Em0 = results(ip).direct.mpcs(5).avg_s_t{1,1} ...
                        - results(baseind).direct.mpcs(5).avg_s_t{1,1};
        decomp2(ip).mpc1_term1 = g0_a' * (m1_a - m0_a); % Effect of MPC function
        decomp2(ip).mpc1_term2 = m0_a' * (g1_a - g0_a); % Effect of distribution
        decomp2(ip).mpc1_term3 = (m1_a - m0_a)' * (g1_a - g0_a); % Interaction

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
	            decomp2(ip).mpc1_term2a(ia) = m0_a(idx)' * (g1_a(idx) - g0_a(idx));
	            % Non-HtM households
	            decomp2(ip).mpc1_term2b(ia) = m0_a(~idx)' * (g1_a(~idx) - g0_a(~idx));
	        else
	        	decomp2(ip).mpc1_term2a(ia) = m0g1_interp(abar) - m0g0_interp(abar);
	        	decomp2(ip).mpc1_term2b(ia) = (m0_a'*g1_a - m0g1_interp(abar)) ...
	        									- (m0_a'*g0_a - m0g0_interp(abar));
	        end
        end
        
        %% 1-Period MPC decomposition3 (decomp with respect to representative agent model (RA))
        if params(ip).nb == 1
            % RA MPC
            m0 = params(ip).R * (results(ip).direct.beta*params(ip).R)^(-1/params(ip).risk_aver) - 1;

            % To get E[MPC|a=3.5], interpolate
            m1_a_interp = griddedInterpolant(agrid_short,m1_a,'linear');
            m1_atmean = m1_a_interp(3.5);

            % Mean model MPC - RA MPC
            decomp3(ip).mpc1_Em1_less_mRA = results(ip).direct.mpcs(5).avg_s_t{1,1} - m0;

            decomp3(ip).mpc1_term1 = m1_atmean - m0; % Effect of MPC function
            decomp3(ip).mpc1_term2 = 0; % Effect of distribution
            decomp3(ip).mpc1_term3 = (m1_a - m0)' * g1_a - (m1_atmean - m0); % Interaction
        end
        

        %% 4-Period MPC decomp2 (decomp with respect to baseline)
        % Only used for quarterly models
        if params(ip).freq == 4
            m1_4 = 0;
            for i = 1:4
                % model MPCs
                m1_it = results(ip).direct.mpcs(5).mpcs_1_t{i};
                m1_4 = m1_4 + m1_it;
            end
            
            m1_4_wide = reshape(m1_4,[params(ip).nx_KFE dim1]);
            m1_4_a = m1_4_wide .* g1_wide ./ g1_a_nonzero;
            m1_4_a = sum(m1_4_a,2);
            
            % P(a) too small to divide by accurately, use arithmetic mean
            m1_4_a(g1small) = mean(m1_4_wide(g1small,:),2);

            % create interpolants from assets to m0g0, m0g1
	        m0g0_4 = m0_4_a .* g0_a;
	        m0g0_interp4 = griddedInterpolant(agrid_short,cumsum(m0g0_4),'linear');
	        m0g1_4 = m0_4_a .* g1_a;
	        m0g1_interp4 = griddedInterpolant(agrid_short,cumsum(m0g1_4),'linear');

            decomp2(ip).mpc4_Em1_less_Em0 = results(ip).direct.mpcs(5).avg_1_1to4 ...
                        - results(baseind).direct.mpcs(5).avg_1_1to4;
            decomp2(ip).mpc4_term1 = g0_a' * (m1_4_a - m0_4_a);
            decomp2(ip).mpc4_term2 = m0_4_a' * (g1_a - g0_a);
            decomp2(ip).mpc4_term3 = (m1_4_a - m0_4_a)' * (g1_a - g0_a);

            for ia = 1:numel(params(ip).abars)
                abar = params(ip).abars(ia);
                if abar == 0
	                idx = agrid_short <= abar;
	                decomp2(ip).mpc4_term2a(ia) = m0_4_a(idx)' * (g1_a(idx) - g0_a(idx));
	                decomp2(ip).mpc4_term2b(ia) = m0_4_a(~idx)' * (g1_a(~idx) - g0_a(~idx));
	            else
	            	decomp2(ip).mpc4_term2a(ia) = m0g1_interp4(abar) - m0g0_interp4(abar);
		        	decomp2(ip).mpc4_term2b(ia) = (m0_4_a'*g1_a - m0g1_interp4(abar)) ...
	        									- (m0_4_a'*g0_a - m0g0_interp4(abar))
            	end
            end
        end
            

    end
end