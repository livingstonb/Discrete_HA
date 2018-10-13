function [distribution,sav] = find_stationary_adist(p,model,...
                                 	income,prefs,agridinput)
    % Finds the stationary distribution and transition matrix for a given
    % agridinput

    if p.Display == 1
        fprintf(' Computing state-to-state transition probabilities... \n');
    end

    gridsize = size(agridinput,1);
 
    NN = gridsize * p.nyP * p.nyF * p.nb;
    nn = gridsize;
    netymat_fulldim = reshape(income.netymat,[p.nyP p.nyF p.nyT]);
    
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
        statetrans_live = bsxfun(@times,transcol_live,interp_live);
        transcol_death = kron(trans_death(:,col),ones(gridsize,1));
        statetrans_death = bsxfun(@times,transcol_death,interp_death);
   
        % add new column to transition matrix
        statetrans(:,nn*(col-1)+1:nn*col) = ...
            (1-p.dieprob)*statetrans_live + p.dieprob*statetrans_death;
        col = col + 1;
    end
    end
    end

    % stationary distribution over states
    if p.Display == 1
        fprintf(' Finding ergodic distribution...\n');
    end
    %distribution = full(ergodicdist(statetrans));
    opts.tol = 1e-8;
    [distribution,~] = eigs(statetrans',1,1,opts);
    distribution = full(distribution/sum(distribution));
    
    distribution = reshape(distribution,[nn,p.nyP,p.nyF,p.nb]);
end
