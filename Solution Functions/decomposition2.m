function [decomp2,decomp3] = decomposition2(params,results)
    % decomp2 is relative to baseline
    % decomp3 is relative to RA model
    
    % computed with respect to a shock of 0.01 * mean ann gross income

    % Construct agrid based off p(1)
    agrid = linspace(0,1,params(1).nxlong)';
    agrid = agrid.^(1/params(1).xgrid_par);
    agrid = params(1).borrow_lim + (params(1).xmax - params(1).borrow_lim) * agrid;

    % Force grid spacing >= gridspace_min near 0
    for ia = 1:params(1).nxlong-1
        if agrid(ia+1) - agrid(ia) < params(1).gridspace_min
            agrid(ia+1) = agrid(ia) + params(1).gridspace_min;
        else
            break
        end
    end

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

        decomp2(ip).mpc1_term3a = NaN(numel(params(1).abars),1);
        decomp2(ip).mpc1_term3b = NaN(numel(params(1).abars),1);
        decomp2(ip).mpc4_term3a = NaN(numel(params(1).abars),1);
        decomp2(ip).mpc4_term3b = NaN(numel(params(1).abars),1);

        decomp2(ip).decomp_error = NaN;
        decomp3(ip).decomp_error = NaN;

        try
            % find index of baseline within 'p'
            if params(ip).freq == 1
                baseind = find(ismember({params.name},{'baseline_A'}));
            else
                baseind = find(ismember({params.name},{'baseline_Q'}));
            end
            
            g1 = results(ip).direct.adist(:);       % model distribution
            g0 = results(baseind).direct.adist(:);  % baseline distribution

            %% 1-Period MPC decomposition2 (decomp with respect to baseline)
            m1 = results(ip).direct.mpcs(5).mpcs_1_t{1};        % model MPCs
            m0 = results(baseind).direct.mpcs(5).mpcs_1_t{1};   % baseline MPCs

            % Mean model MPC  - mean baseline MPC
            decomp2(ip).mpc1_Em1_less_Em0 = results(ip).direct.mpcs(5).avg_s_t{1,1} ...
                            - results(baseind).direct.mpcs(5).avg_s_t{1,1};
            decomp2(ip).mpc1_term1 = g0' * (m1 - m0); % Effect of MPC function
            decomp2(ip).mpc1_term2 = m0' * (g1 - g0); % Effect of distribution
            decomp2(ip).mpc1_term3 = (m1 - m0)' * (g1 - g0); % Interaction

            % Decomposition of distribution effect
            for ia = 1:numel(params(ip).abars)
                abar = params(ip).abars(ia); % HtM threshold
                idx = agrid <= abar;
                % HtM households
                decomp2(ip).mpc1_term3a(ia) = m0(idx)' * (g1(idx) - g0(idx));
                % Non-HtM households
                decomp2(ip).mpc1_term3b(ia) = m0(~idx)' * (g1(~idx) - g0(~idx));
            end
            
            %% 1-Period MPC decomposition3 (decomp with respect to representative agent model (RA))
            if params(ip).nb == 1
                % model MPCs
                m1 = results(ip).direct.mpcs(5).mpcs_1_t{1};
                % RA MPC
                m0 = params(ip).R * (results(ip).direct.beta*params(ip).R)^(-1/params(ip).risk_aver) - 1;

                % expected value of m1 conditional on assets = 0
                m1_at0 = m1(agrid(:)==0)' * g1(agrid(:)==0) / sum(g1(agrid(:)==0));

                % Mean model MPC - RA MPC
                decomp3(ip).mpc1_Em1_less_mRA = results(ip).direct.mpcs(5).avg_s_t{1,1} - m0;

                decomp3(ip).mpc1_term1 = m1_at0- m0; % Effect of MPC function
                decomp3(ip).mpc1_term2 = 0; % Effect of distribution
                decomp3(ip).mpc1_term3 = (m1 - m0)' * g1 - m1_at0 - m0; % Interaction
            end
            

            %% 4-Period MPC decomp2 (decomp with respect to baseline)
            % Only used for quarterly models
            if params(ip).freq == 4
                m1 = 0;
                m0 = 0;
                for i = 1:4
                    % model MPCs
                    m1 = m1 + results(ip).direct.mpcs(5).mpcs_1_t{i};
                    % baseline MPCs
                    m0 = m0 + results(baseind).direct.mpcs(5).mpcs_1_t{i};
                end

                decomp2(ip).mpc4_Em1_less_Em0 = results(ip).direct.mpcs(5).avg_1_1to4 ...
                            - results(baseind).direct.mpcs(5).avg_1_1to4;
                decomp2(ip).mpc4_term1 = g0' * (m1 - m0);
                decomp2(ip).mpc4_term2 = m0' * (g1 - g0);
                decomp2(ip).mpc4_term3 = (m1 - m0)' * (g1 - g0);

                for ia = 1:numel(params(ip).abars)
                    abar = params(ip).abars(ia);
                    idx = agrid <= abar;
                    decomp2(ip).mpc4_term3a(ia) = m0(idx)' * (g1(idx) - g0(idx));
                    decomp2(ip).mpc4_term3b(ia) = m0(~idx)' * (g1(~idx) - g0(~idx));
                end
            end
            
            
        catch ME
            decomp2(ip).decomp_error = ME;
            decomp3(ip).decomp_error = ME;
        end
    end