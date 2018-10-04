function [distribution,statetrans,sav,con] = find_stationary(p,model,...
                            income,prefs,xgridinput,ergodic_method,ergodic_tol)
    % Finds the stationary distribution and transition matrix for a given
    % xgridinput. sav and con are the saving and consumption policy functions
    % associated with xgridinput
    
    % xgridinput should be of the _wide format: dimensions (assets,nyP,nyF)
    
    % ergodic_method is 1 for iterative, 2 for direct

    if p.Display == 1
        fprintf(' Computing state-to-state transition probabilities... \n');
    end

    gridsize = size(xgridinput,1);
    
    % Interpolate policy functions onto xgridinput
    sav = zeros(gridsize,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        sav(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(xgridinput(:,iyP,iyF));
        con(:,iyP,iyF,ib) = model.coninterp{iyP,iyF,ib}(xgridinput(:,iyP,iyF));
    end
    end
    end
    
    NN = gridsize * p.nyP * p.nyF * p.nb;
    nn = gridsize;
    netymat_fulldim = reshape(income.netymat,[p.nyP p.nyF p.nyT]);
    
    % transition matrix between (yP,yF,beta) states, condl on living
    trans_live = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    % find transtion matrix between (yP,yF,beta) states condl on dying
    ytrans_stationary = repmat(income.yPdist',p.nyP,1);
    trans_death = kron(prefs.betatrans,kron(eye(p.nyF),ytrans_stationary));
    
    % transition matrix over (x,yP,yF,beta) full asset space
    statetrans = zeros(NN,NN);
    col = 1;
    for ib2 = 1:p.nb
    for iyF2 = 1:p.nyF
    for iyP2 = 1:p.nyP
        % create spline object
        fspace = fundef({'spli',xgridinput(:,iyP2,iyF2),0,1});
        
        % xprime if no death
        xp = (1+p.r)*repmat(sav(:),p.nyT,1) + ...
            kron(squeeze(netymat_fulldim(iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
        
        xp_live  = xp;
        if p.WealthInherited == 0
            % x' is only equal to income
            xp_death = kron(squeeze(netymat_fulldim(iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
        else
            xp_death = xp_live;
        end
        
        % xp_live and xp_death interpolated onto xgridinput
        interp_live     = reshape(funbas(fspace,xp_live),[],p.nyT*nn);
        interp_death    = reshape(funbas(fspace,xp_death),[],p.nyT*nn);

        % take expectation over yT distribution
        interp_live     = interp_live  * kron(speye(nn),income.yTdist);
        interp_death    = interp_death * kron(speye(nn),income.yTdist);

        % Multiply by transition matrix between (yP,yF,beta) states
        interp_live     = bsxfun(@times,kron(trans_live(:,col),ones(nn,1)),interp_live);
        interp_death    = bsxfun(@times,kron(trans_death(:,col),ones(nn,1)),interp_death);

        % add new column to transition matrix
        newcolumn       = (1-p.dieprob)*interp_live + p.dieprob*interp_death;
        statetrans(:,nn*(col-1)+1:nn*col) = newcolumn;
        col = col + 1;
    end
    end
    end
    
    statetrans = sparse(statetrans);

    % stationary distribution over states
    if ergodic_method == 0
        distribution = [];
    else
        fprintf(' Finding ergodic distribution...\n');
        distribution = double(full(ergodicdist(statetrans,ergodic_method,ergodic_tol)));
    end
end
