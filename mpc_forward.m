function [avg_mpc1,avg_mpc4,gridspace] = mpc_forward(xgrid,p,income,model,...
                                                            prefs)
    ergodic_tol = 1e-7;
    
    disp('Using direct methods to find MPCs')
    
    %% CREATE NEW GRID
    nn = p.nxmpc;
    NN = nn * p.nyP * p.nyF * p.nb;
    % Put 1/2 of grid pts in bottom 5% of grid
    xgridm = linspace(0,1,p.nxmpc)';
    xgridm = p.borrow_lim + (p.xmax-p.borrow_lim)*xgridm;
    gridspace = xgridm(2) - xgridm(1);

    % get interpolated saving function on grid
    savm = zeros(nn,p.nyP,p.nyF,p.nb);
    savinterp = cell(p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        savm(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(xgridm);
    end
    end
    end

    % create transition probability matrix
    trans = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    T = transition(income,NN,nn,p,savm,xgridm,trans);

    % find stationary distribution
    stat_dist = full(ergodicdist(sparse(T),1,ergodic_tol));

    %% SHIFT DISTRIBUTION (GIVE CASH TO EVERYONE)
    nn = nn + 1;
    NN = nn * p.nyP * p.nyF * p.nb;
    stat_dist_wide = reshape(stat_dist,[nn-1 p.nyP p.nyF p.nb]);
    stat_dist = zeros(nn,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        stat_dist(2:end,iyP,iyF,ib) = stat_dist_wide(:,iyP,iyF,ib);
        stat_dist(1,iyP,iyF,ib) = 0;
    end
    end
    end
    stat_dist = stat_dist(:);
    xgridm(end+1) = xgridm(end) + gridspace;

    % get interpolated saving function on new grid
    savm = zeros(nn,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        savm(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(xgridm);
    end
    end
    end

    % find new transition probability matrix
    T = transition(income,NN,nn,p,savm,xgridm,trans);

    % find new stationary distribution
    SS_dist = full(ergodicdist(T,1,ergodic_tol));

    % find future average consumption
    savm = zeros(nn,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        savm(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(xgridm);
    end
    end
    end
    savm(savm<p.borrow_lim) = p.borrow_lim;
    conm = repmat(xgridm,[1 p.nyP p.nyF p.nb]) - savm - p.savtax * max(savm - p.savtaxthresh,0);

    %% FIND MPCs

    dist{1} = stat_dist';
    dist{2} = dist{1} * T;
    dist{3} = dist{2} * T;
    dist{4} = dist{3} * T;

    % Later, use multiples of gridspace
    mpcamount = gridspace;
    avg_mpc{1} = (dist{1} - SS_dist') * conm(:)/mpcamount;

    for period = 2:4
        avg_mpc{period} = avg_mpc{period-1} + (dist{period} - SS_dist') * conm(:)/mpcamount;
    end
    
    avg_mpc1 = avg_mpc{1};
    avg_mpc4 = avg_mpc{4};

    %% FUNCTIONS
    
    function T = transition(income,NN,nn,p,savm,xgridm,trans)
        fspace = fundef({'spli',xgridm,0,1});
        % xprime conditional on dying
        xp_death = repmat(income.netymat(:)',NN,1);
        % xprime conditional on living
        xp_live = xp_death + p.R * repmat(savm(:),1,p.nyT*p.nyF*p.nyP);

        % interpolated conditional transition probabilities
        interp_death = funbas(fspace,xp_death(:));
        interp_live = funbas(fspace,xp_live(:));
        % interpolated transition probabilities
        interp = (1-p.dieprob)*interp_live + p.dieprob*interp_death;
        % take expected value over yT
        interp = reshape(interp,[],nn*p.nyT) * kron(speye(nn),income.yTdist);
        % reshape and multiply by transition probabilities of (yP,yF,beta)
        interp = permute(reshape(full(interp),[NN p.nyP*p.nyF nn]),[1 3 2]);
        interp = reshape(interp,NN,nn*p.nyP*p.nyF);
        T = sparse(kron(trans,ones(nn)) .* repmat(interp,1,p.nb));
    end

end