classdef DecompTwo < handle
    
    properties
        mpc1_Em1_less_Em0;
        mpc1_term1;
        mpc1_term2;
        mpc1_term3;
        mpc1_term3a;
        mpc1_term3b;
        
        mpc4_Em1_less_Em0;
        mpc4_term1;
        mpc4_term2;
        mpc4_term3;
        mpc4_term3a;
        mpc4_term3b;
    end

    methods
        function obj = DecompTwo(params)
            % Initialize all fields to NaN
            fields = fieldnames(obj);
            for i = 1:numel(fields)
                obj.(fields{i}) = NaN;
            end
            
            % Set baseline index
            
        end
        
    end
    
    methods (Static)
        function objs = decompose(objs,params,results)
            % Construct agrid based off params(1)
            agrid = linspace(0,1,params(1).nxlong)';
            agrid = agrid.^(1/params(1).xgrid_par);
            agrid = params(1).borrow_lim + (params(1).xmax - params(1).borrow_lim) * agrid;
            
            for ip = 1:numel(objs)
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
            end  
        end
    end
    
end