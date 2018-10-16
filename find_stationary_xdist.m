function [distribution,sav] = find_stationary_xdist(p,model,...
                                 	income,prefs,xgridinput)
    % Finds the stationary distribution and transition matrix for a given
    % xgridinput. sav is the saving policy function
    % associated with xgridinput.

    if p.Display == 1
        fprintf(' Computing state-to-state transition probabilities... \n');
    end

    gridsize = size(xgridinput,1);
    
    if ismatrix(xgridinput)==1 && size(xgridinput,2)==1
        xgridinput = repmat(xgridinput,[1 p.nyP p.nyF]);
    end
    
    % Interpolate policy functions onto xgridinput
    sav = zeros(gridsize,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        sav(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(xgridinput(:,iyP,iyF));
    end
    end
    end
    
    NN = gridsize * p.nyP * p.nyF * p.nb;
    nn = gridsize;
    netymat_fulldim = reshape(income.netymat,[p.nyP p.nyF p.nyT]);
    
    % transition matrix between (yP,yF,beta) states, cond'l on living
    trans_live = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    % transition matrix between (yP,yF,beta) states cond'l on dying
    yPtrans_stationary = repmat(income.yPdist',p.nyP,1);
    trans_death = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans_stationary));
    
    % transition matrix over (x,yP,yF,beta) full asset space
    statetrans = sparse(NN,NN);
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
        if p.Bequests == 0
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
        trans_live_expand = kron(trans_live(:,col),ones(nn,1));
        interp_live     = bsxfun(@times,trans_live_expand,interp_live);
        trans_death_expand = kron(trans_death(:,col),ones(nn,1));
        interp_death    = bsxfun(@times,trans_death_expand,interp_death);

        % add new column to transition matrix
        statetrans(:,nn*(col-1)+1:nn*col) = (1-p.dieprob)*interp_live...
                                                + p.dieprob*interp_death;
        col = col + 1;
    end
    end
    end
    clear interp_live interp_death trans_live_expand trans_death_expand

    % stationary distribution over states
    if p.Display == 1
        fprintf(' Finding ergodic distribution...\n');
    end

    if p.nyF == 1 && (p.nb==1 ||( p.nb>1 && p.betaswitch>0))
        % No fixed heterogeneity
        opts.v0 = sparse(NN,1);
        opts.v0(1) = 1;
        [distribution,~] = eigs(statetrans',1,1,opts);
        distribution = distribution/sum(distribution);
    else
        % Need iterative procedure
        q=sparse(1,NN);
        % Create valid initial distribution for both yF & beta
        % Repmat automatically puts equal weight on each beta
        q(1,1:nn*p.nyP:end)=repmat(income.yFdist,p.nb,1);
        diff=1; 
        iter = 1;
        while diff>1e-8 && iter < 1e6;
            z=q*statetrans;
            diff=norm(z-q);
            q=z;
            iter = iter + 1;
        end
        if iter >= 1e6
            error('No convergence to stationary distribution')
        end
        distribution=full(q');
    end
    
    distribution = reshape(distribution,[nn p.nyP p.nyF p.nb]);

end
