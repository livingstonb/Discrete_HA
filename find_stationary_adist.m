function [adist,xdist,xvals,incvals,netincvals] = find_stationary_adist(p,model,...
                                 	income,prefs,agridinput)
    % Finds the stationary distribution and transition matrix for a given
    % agridinput

    if p.Display == 1
        fprintf(' Computing state-to-state transition probabilities... \n');
    end

    gridsize = size(agridinput,1);
 
    NN = gridsize * p.nyP * p.nyF * p.nb;
    nn = gridsize;
    
    % transition matrix between (yP,yF,beta) states, cond'l on living
    trans_live = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    % transition matrix between (yP,yF,beta) states cond'l on dying
    yPtrans_stationary = repmat(income.yPdist',p.nyP,1);
    trans_death = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans_stationary));
    
    agrid_full = repmat(agridinput,[1 p.nyP p.nyF p.nyT]);
    netymat_full = reshape(income.netymat,[1 p.nyP p.nyF p.nyT]);
    netymat_full = repmat(netymat_full,[p.nxlong 1 1 1]);
    
    % cash-on-hand as function of (a,yP,yF,yT)
    x = agrid_full + netymat_full;
    
    % saving interpolated onto this grid
    sav = zeros(gridsize,p.nyP,p.nyF,p.nb,p.nyT);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        x_iyP_iyF = x(:,iyP,iyF,:);
        sav_iyP_iyF_ib = model.savinterp{iyP,iyF,ib}(x_iyP_iyF(:));
        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_ib,[p.nxlong 1 1 1 p.nyT]);
    end
    end
    end
    
    aprime_live = p.R * sav;

    % transition matrix over (x,yP,yF,beta) full asset space
    statetrans = sparse(NN,NN);
    % create spline object
    fspace = fundef({'spli',agridinput,0,1});
    % get interpolated probabilities and take expectation over yT
    interp_live = funbas(fspace,aprime_live(:));
    interp_live = reshape(interp_live,NN,nn*p.nyT);
    interp_live = interp_live * kron(speye(nn),income.yTdist);
    if p.Bequests == 1
        interp_death = interp_live;
    else
        interp_death = sparse(NN,nn);
        interp_death(:,1) = ones(NN,1);
    end

    col = 1;
    for ib2 = 1:p.nb
    for iyF2 = 1:p.nyF
    for iyP2 = 1:p.nyP
        transcol_live = kron(trans_live(:,col),ones(gridsize,1));
        transcol_live = bsxfun(@times,transcol_live,interp_live);
        transcol_death = kron(trans_death(:,col),ones(gridsize,1));
        transcol_death = bsxfun(@times,transcol_death,interp_death);
   
        % add new column to transition matrix
        statetrans(:,nn*(col-1)+1:nn*col) = ...
            (1-p.dieprob)*transcol_live + p.dieprob*transcol_death;
        col = col + 1;
    end
    end
    end
    clear transcol_live transcol_death interp_live interp_death

    % stationary distribution over states
    if p.Display == 1
        fprintf(' Finding ergodic distribution...\n');
    end
    
    if p.nyF == 1 && (p.nb==1 ||( p.nb>1 && p.betaswitch>0))
        % No fixed heterogeneity
        opts.v0 = sparse(NN,1);
        opts.v0(1) = 1;
        [adist,~] = eigs(statetrans',1,1,opts);
        adist = adist/sum(adist);
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
        adist=full(q');
    end
    
    adist = reshape(adist,[nn,p.nyP,p.nyF,p.nb]);
    
    % get distribution over (x,yP,yF,beta)
    xdist = kron(income.yTdist,reshape(adist,nn,[]));
    xdist = reshape(xdist,[nn*p.nyT p.nyP p.nyF p.nb]);
    
    % Extend xvals to (p.nxlong*p.nyT,p.nyP,p.nyF,p.nyT)
    incvals = reshape(income.ymat,[p.nyP*p.nyF p.nyT]);
    incvals = permute(incvals,[2 1]);
    incvals = kron(incvals,ones(nn,1));
    incvals = reshape(incvals,[nn*p.nyT p.nyP p.nyF]);
    incvals = repmat(incvals,[1 1 1 p.nb]);
    netincvals = income.lumptransfer + (1-p.labtaxlow)*incvals - p.labtaxhigh*max(incvals-income.labtaxthresh,0);
    xvals = repmat(agridinput,[p.nyT p.nyP p.nyF p.nb]) + netincvals;
end
