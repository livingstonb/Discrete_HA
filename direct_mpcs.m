function [mpcs1,mpcs4,avg_mpc1,avg_mpc4] = direct_mpcs(xgrid,p,income,basemodel,prefs)
    % Uses direct methods to compute MPCs
    
    if p.Display == 1
        disp('Using direct methods to compute MPCs...')
    end
    
    NN          = p.nxlong * p.nyP * p.nyF * p.nb;
    mpcs1      = cell(1,numel(p.mpcfrac));
    mpcs4      = cell(1,numel(p.mpcfrac));
    avg_mpc1   = cell(1,numel(p.mpcfrac));
    avg_mpc4   = cell(1,numel(p.mpcfrac));
    
    statetrans2 = basemodel.statetrans^2;
    
    for im = 1:numel(p.mpcfrac)
        mpcamount = p.mpcfrac{im} * income.meany * p.freq;
        x_mpc = xgrid.longgrid_wide + mpcamount;
        
        % Get policy functions on shifted grid
        con = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        sav = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            con(:,iyP,iyF,ib) = basemodel.coninterp{iyP,iyF,ib}(x_mpc(:,iyP,iyF));
            sav(:,iyP,iyF,ib) = basemodel.savinterp{iyP,iyF,ib}(x_mpc(:,iyP,iyF));
        end
        end
        end
        
        % One-period MPC
        mpcs{1} = (con(:) - basemodel.con_longgrid) / mpcamount;
        avg_mpc{1} = basemodel.SSdist' * mpcs{1};

        % Create matrix of xprime's, final dim NN by p.nyP*p.nyF*p.nyT
        % Rows index points in today's asset space, columns index
        % next period's yP,yF,yT realizations
        xprime   = zeros(NN,p.nyP*p.nyF,p.nyT);
        for iyT = 1:p.nyT
            xprime(:,:,iyT) = p.R * repmat(sav(:),1,p.nyP*p.nyF)...
                                    + repmat(income.netymat(:,iyT)',NN,1);
        end
        xprime = reshape(xprime,NN,p.nyP*p.nyF*p.nyT);

        % Get transition matrix
        T = transition_matrix(NN,xprime,p,xgrid,income,prefs);

        % Multi-period MPCs
        mpcs{2} = (T * basemodel.con_longgrid - basemodel.con_longgrid) / mpcamount;
        mpcs{3} = (T * basemodel.statetrans  * basemodel.con_longgrid - basemodel.con_longgrid) / mpcamount;
        mpcs{4} = (T * statetrans2 * basemodel.con_longgrid - basemodel.con_longgrid) / mpcamount;
        
        for it = 2:4;
            avg_mpc{it} = avg_mpc{it-1} + basemodel.SSdist' * mpcs{it};
        end
        
        % Store one- and four-period MPCs
        mpcs1{im}       = mpcs{1};
        mpcs4{im}       = mpcs{4};
        avg_mpc1{im}    = avg_mpc{1};
        avg_mpc4{im}    = avg_mpc{4};
    end   

    function T = transition_matrix(NN,xp,p,xgrid,income,prefs)
        T = sparse(NN,NN);
        col = 1;
        xp = reshape(xp,NN,p.nyP,p.nyF,p.nyT);
        for ib2 = 1:p.nb
        for iyF2 = 1:p.nyF
        for iyP2 = 1:p.nyP
            % Create spline object
            fspace = fundef({'spli',xgrid.longgrid_wide(:,iyP2,iyF2),0,1});

            % Interpolate x' onto the grid
            xp_T = reshape(squeeze(xp(:,iyP2,iyF2,:)),[],1);
            newcolumn = funbas(fspace,xp_T);
            
            % Take expectation over yT distribution
            newcolumn = reshape(newcolumn,[],p.nyT*p.nxlong) * kron(speye(p.nxlong),income.yTdist);
            
            T(:,p.nxlong*(col-1)+1:p.nxlong*col) = newcolumn;
            col = col + 1;
        end
        end
        end
        
        % Transition matrix of exogenous states
        trans_exo = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
        trans_exo = kron(trans_exo,ones(p.nxlong));
        
        T = trans_exo .* T;
    end

end