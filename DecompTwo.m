classdef DecompTwo < handle
    
    properties
        mpc1_Em1_less_Em0 = NaN;
        mpc1_term1 = NaN;
        mpc1_term2 = NaN;
        mpc1_term3 = NaN;
        mpc1_term3a;
        mpc1_term3b;
        
        mpc4_Em1_less_Em0 = NaN;
        mpc4_term1 = NaN;
        mpc4_term2 = NaN;
        mpc4_term3 = NaN;
        mpc4_term3a;
        mpc4_term3b;
        
        decomp_error;
    end

    methods
        function obj = DecompTwo(params_ip)
            
            obj.mpc1_term3a = NaN(numel(params_ip.abars),1);
            obj.mpc1_term3b = NaN(numel(params_ip.abars),1);
            obj.mpc4_term3a = NaN(numel(params_ip.abars),1);
            obj.mpc4_term3b = NaN(numel(params_ip.abars),1);
        end
        
    end
    
    methods (Static)
        function objs = decompose(objs,params,results)
            % Construct agrid based off params(1)
            agrid = linspace(0,1,params(1).nxlong)';
            agrid = agrid.^(1/params(1).xgrid_par);
            agrid = params(1).borrow_lim + (params(1).xmax - params(1).borrow_lim) * agrid;
            
            for ip = 1:numel(objs)
                try
                    if params(ip).freq == 1
                        baseind = find(ismember({params.name},{'baseline_A'}));
                    else
                        baseind = find(ismember({params.name},{'baseline_Q'}));
                    end

                    g1 = results(ip).direct.agrid_dist;
                    g0 = results(baseind).direct.agrid_dist;

                    % 1-Period MPC decomp
                    m1 = results(ip).direct.mpcs1_a_direct{5};
                    m0 = results(baseind).direct.mpcs1_a_direct{5};

                    objs(ip).mpc1_Em1_less_Em0 = results(ip).direct.avg_mpc1_agrid(5) ...
                                    - results(baseind).direct.avg_mpc1_agrid(5);
                    objs(ip).mpc1_term1 = g0' * (m1 - m0);
                    objs(ip).mpc1_term2 = m0' * (g1 - g0);
                    objs(ip).mpc1_term3 = (m1 - m0)' * (g1 - g0);

                    for ia = 1:numel(params(ip).abars)
                        abar = params(ip).abars(ia);
                        idx = agrid <= abar;
                        objs(ip).mpc1_term3a(ia) = m0(idx)' * (g1(idx) - g0(idx));
                        objs(ip).mpc1_term3b(ia) = m0(~idx)' * (g1(~idx) - g0(~idx));
                    end

                    % 4-Period MPC decomp
                    if params(ip).freq == 4
                        m1 = results(ip).direct.mpcs4_a_direct{5};
                        m0 = results(baseind).direct.mpcs4_a_direct{5};

                        objs(ip).mpc4_Em1_less_Em0 = results(ip).direct.avg_mpc4_agrid(5) ...
                                    - results(baseind).direct.avg_mpc4_agrid(5);
                        objs(ip).mpc4_term1 = g0' * (m1 - m0);
                        objs(ip).mpc4_term2 = m0' * (g1 - g0);
                        objs(ip).mpc4_term3 = (m1 - m0)' * (g1 - g0);

                        for ia = 1:numel(params(ip).abars)
                            abar = params(ip).abars(ia);
                            idx = agrid <= abar;
                            objs(ip).mpc4_term3a(ia) = m0(idx)' * (g1(idx) - g0(idx));
                            objs(ip).mpc4_term3b(ia) = m0(~idx)' * (g1(~idx) - g0(~idx));
                        end
                    end
                catch ME
                    objs(ip).decomp_error = ME;
                end
            end  
        end
    end
    
end