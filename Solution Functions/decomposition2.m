function [decomp2,decomp3] = decomposition2(params,results)
    % decomp2 is relative to baseline
    % decomp3 is relative to RA model

    % Construct agrid based off p(1)
    agrid = linspace(0,1,params(1).nxlong)';
    agrid = agrid.^(1/params(1).xgrid_par);
    agrid = params(1).borrow_lim + (params(1).xmax - params(1).borrow_lim) * agrid;
    % Force grid spacing >= gridspace_min near 0, assumes all gridspace_min
    % equal
    for ia = 1:params(1).nxlong-1
        if agrid(ia+1) - agrid(ia) < params(1).gridspace_min
            agrid(ia+1) = agrid(ia) + params(1).gridspace_min;
        else
            break
        end
    end
    
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
            if params(ip).freq == 1
                baseind = find(ismember({params.name},{'baseline_A'}));
            else
                baseind = find(ismember({params.name},{'baseline_Q'}));
            end
            
            g1 = results(ip).direct.adist(:);
            g0 = results(baseind).direct.adist(:);

            % 1-Period MPC decomp2
            m1 = results(ip).direct.mpcs01.mpcs_1_t{1};
            m0 = results(baseind).direct.mpcs01.mpcs_1_t{1};

            decomp2(ip).mpc1_Em1_less_Em0 = results(ip).direct.mpcs01.avg_s_t{1,1} ...
                            - results(baseind).direct.mpcs01.avg_s_t{1,1};
            decomp2(ip).mpc1_term1 = g0' * (m1 - m0);
            decomp2(ip).mpc1_term2 = m0' * (g1 - g0);
            decomp2(ip).mpc1_term3 = (m1 - m0)' * (g1 - g0);

            for ia = 1:numel(params(ip).abars)
                abar = params(ip).abars(ia);
                idx = agrid <= abar;
                decomp2(ip).mpc1_term3a(ia) = m0(idx)' * (g1(idx) - g0(idx));
                decomp2(ip).mpc1_term3b(ia) = m0(~idx)' * (g1(~idx) - g0(~idx));
            end
            
            % 1-Period MPC decomp3
            if params(ip).nb == 1
                m1 = results(ip).direct.mpcs01.mpcs_1_t{1};
                m0 = params(ip).R * (results(ip).direct.beta*params(ip).R)^(-1/params(ip).risk_aver) - 1;


                decomp3(ip).mpc1_Em1_less_mRA = results(ip).direct.mpcs01.avg_s_t{1,1} - m0;
                decomp3(ip).mpc1_term1 = m1(1)- m0; % m1 at assets = 0
                decomp3(ip).mpc1_term2 = 0;
                decomp3(ip).mpc1_term3 = (m1 - m0)' * g1 - m1(1) - m0;
            end
            

            % 4-Period MPC decomp
            if params(ip).freq == 4
                m1 = 0;
                m0 = 0;
                for i = 1:4
                    m1 = m1 + results(ip).direct.mpcs01.mpcs_1_t{i};
                    m0 = m0 + results(baseind).direct.mpcs01.mpcs_1_t{i};
                end

                decomp2(ip).mpc4_Em1_less_Em0 = results(ip).direct.mpcs01.avg_1_1to4 ...
                            - results(baseind).direct.mpcs01.avg_1_1to4;
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