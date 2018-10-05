function [mpcs1,mpcs4,avg_mpc1,avg_mpc4,norisk] = expected_mpcs_deterministic(...
                                                xgrid,p,income,norisk,prefs);
    % Uses direct methods to compute MPCs
    
    if p.Display == 1
        disp('Using direct methods to compute MPCs (norisk model)...')
    end
    
    NN          = p.nxlong * p.nb;
    mpcs1      = cell(1,numel(p.mpcfrac));
    mpcs4      = cell(1,numel(p.mpcfrac));
    avg_mpc1   = cell(1,numel(p.mpcfrac));
    avg_mpc4   = cell(1,numel(p.mpcfrac));
    
    % interpolate policy functions onto longer grid
    norisk.con_longgrid_wide = zeros(p.nxlong,p.nb);
    norisk.sav_longgrid_wide = zeros(p.nxlong,p.nb);
    for ib = 1:p.nb
        norisk.con_longgrid_wide(:,ib) = norisk.coninterp{ib}(xgrid.norisk_longgrid);
        norisk.sav_longgrid_wide(:,ib) = norisk.coninterp{ib}(xgrid.norisk_longgrid);
    end
    
    % find transition matrix and stationary distribution
    trans_exo_norisk = kron(prefs.betatrans,eye(p.nxlong));
    norisk.SSdist_wide = zeros(p.nxlong,p.nb);
    
    fspace = fundef({'spli',xgrid.norisk_longgrid,0,1});
    
    if p.WealthInherited == 1
        xp = p.R * norisk.sav_longgrid_wide + income.meany;
    else
        xp_life  = p.R * norisk.sav_longgrid_wide + income.meany;
        xp_death = income.meany * ones(p.nxlong,p.nb);
    end
    
    norisk.statetrans = sparse(p.nxlong*p.nb,p.nxlong*p.nb);
    for ib = 1:p.nb
        if p.WealthInherited == 1
            interp = funbas(fspace,xp(:));
        else
            interp_life = funbas(fspace,xp_life(:,ib));
            interp_death = funbas(fspace,xp_death(:,ib));
            interp = (1-p.dieprob)*interp_life + p.dieprpob*interp_death;
        end
        norisk.statetrans(:,p.nxlong*(ib-1)+1:p.nxlong*ib) = interp;
    end
    norisk.statetrans = trans_exo_norisk .* norisk.statetrans;
    norisk.SSdist = full(ergodicdist(norisk.statetrans,1,1e-8));

    % Loop over mpcamounts
    for im = 1:numel(p.mpcfrac)
        mpcamount = p.mpcfrac{im} * income.meany * p.freq;
        x_mpc = xgrid.norisk_longgrid + mpcamount;
        
        % if shock is negative, deal with households that wind up
        % below the bottom of their asset grid
        if mpcamount < 0
            set_mpc_one = false(p.nxlong,p.nb);
            for ib = 1:p.nb
                set_mpc_one(:,ib) = x_mpc(:,ib) < xgrid.longgrid_wide(1,ib);
                x_mpc(x_mpc(:,ib)<xgrid.longgrid_wide(1,ib)) = xgrid.longgrid_wide(1,ib);
            end
        end
        
        % Get policy functions on shifted grid
        con = zeros(p.nxlong,p.nb);
        sav = zeros(p.nxlong,p.nb);
        for ib  = 1:p.nb
            con(:,ib) = norisk.coninterp{ib}(x_mpc);
            sav(:,ib) = norisk.savinterp{ib}(x_mpc);
        end

        % One-period MPC
        mpcs{1} = (con(:) - norisk.con_longgrid_wide(:)) / mpcamount;
        if mpcamount < 0
            mpcs{1}(set_mpc_one(:)) = 1;
        end
        avg_mpc{1} = norisk.SSdist' * mpcs{1};

        if p.freq == 1
            % Don't compute four-period MPCs
            mpcs1{im}       = mpcs{1};
            mpcs4{im}       = [];
            avg_mpc1{im}    = avg_mpc{1};
            avg_mpc4{im}    = [];
            continue
        end
        
        % Create matrix of xprime's, final dim p.nxlong*p.nb
        xprime       = p.R * sav(:) + income.meannety * ones(p.nxlong*p.nb,1);
        xprime_death = income.meannety * ones(p.nxlong*p.nb,1);

        % Transition matrix from period t+k to t+k+1 for k > 0
        trans_exo       = kron(prefs.betatrans,eye(p.nxlong));
        
        
        
               
        % Transition matrix from t to t+1
        if p.WealthInherited == 1
            T = trans_exo .* transition_matrix(xprime,p,xgrid);
        else
            T_life  = trans_exo  .* transition_matrix(xprime,p,xgrid);
            T_death = trans_exo  .* transition_matrix(xprime_death,p,xgrid);
            T = (1 - p.dieprob) * T_life + p.dieprob * T_death;
        end

        % Multi-period MPCs
        mpcs{2} = (T * norisk.con_longgrid_wide(:) - norisk.con_longgrid_wide(:)) / mpcamount;
        mpcs{3} = (T * norisk.statetrans   * norisk.con_longgrid_wide(:) - norisk.con_longgrid_wide(:)) / mpcamount;
        mpcs{4} = (T * norisk.statetrans^2 * norisk.con_longgrid_wide(:) - norisk.con_longgrid_wide(:)) / mpcamount;
        
        for it = 2:4;
            avg_mpc{it} = avg_mpc{it-1} + norisk.SSdist' * mpcs{it};
        end
        
        % Store one- and four-period MPCs
        mpcs1{im}       = mpcs{1};
        mpcs4{im}       = mpcs{4};
        avg_mpc1{im}    = avg_mpc{1};
        avg_mpc4{im}    = avg_mpc{4};
    end   

    function T = transition_matrix(xp,p,xgrid)
        T = sparse(p.nxlong*p.nb,p.nxlong*p.nb);
        col = 1;
        for ib2 = 1:p.nb
            % Create spline object
            fspace = fundef({'spli',xgrid.norisk_longgrid,0,1});

            % Interpolate x' onto the grid
            T(:,p.nxlong*(col-1)+1:p.nxlong*col) = funbas(fspace,xp(:));
            col = col + 1;
        end
 
    end

end