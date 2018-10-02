function [distribution,statetrans,sav,con] = find_stationary(p,model,income,prefs,xgridinput,ergodic_tol)
    if p.Display == 1
        fprintf(' Computing state-to-state transition probabilities... \n');
    end
    
    % xgridinput should be of the _wide format: dimensions (assets,nyP,nyF)
    
    
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
    
    trans = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    statetrans = zeros(NN,NN);
    col = 1;
    for ib2 = 1:p.nb
    for iyF2 = 1:p.nyF
    for iyP2 = 1:p.nyP
        fspace = fundef({'spli',xgridinput(:,iyP2,iyF2),0,1});
        % xprime if no death
        xp_live = (1+p.r)*repmat(sav(:),p.nyT,1) + ...
            kron(squeeze(netymat_fulldim(iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));

        % xprime if death (i.e. saving in past period set to 0)
        if p.WealthInherited == 0
            xp_death = kron(squeeze(netymat_fulldim(iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
        else
            xp_death = xp_live;
        end
        
        % set probabilities equal to 1 at grid endpt when xp is off the grid
%         idx_xpl_max = xp_live>=max(xgridinput(:,iyP2,iyF2));
%         idx_xpl_min = xp_live<=min(xgridinput(:,iyP2,iyF2));
%         
%         idx_xpd_max = xp_death>=max(xgridinput(:,iyP2,iyF2));
%         idx_xpd_min = xp_death<=min(xgridinput(:,iyP2,iyF2));
%         
%         interpl = funbas(fspace,xp_live);
%         interpl(idx_xpl_max | idx_xpl_min,:) = 0;
%         interpl(idx_xpl_max,end) = 1;
%         interpl(idx_xpl_min,1) = 1;
%         
%         interpd = funbas(fspace,xp_death);
%         interpd(idx_xpd_max | idx_xpd_min,:) = 0;
%         interpd(idx_xpd_max,end) = 1;
%         interpd(idx_xpd_min,1) = 1;
%         
%         interp = (1-p.dieprob) * interpl + p.dieprob * interpd;
        
            % if not setting probabilites to 1 at grid endpts
        interp = (1-p.dieprob) * funbas(fspace,xp_live) + p.dieprob * funbas(fspace,xp_death);
            
        interp = reshape(interp,[],p.nyT*nn);
        % Multiply by yT distribution
        newcolumn = interp * kron(speye(nn),income.yTdist);
        % Multiply by (beta,yF,yP) distribution
        newcolumn = bsxfun(@times,kron(trans(:,col),ones(nn,1)),newcolumn);

        statetrans(:,nn*(col-1)+1:nn*col) = newcolumn;
        col = col + 1;
    end
    end
    end

    % SS probability of residing in each state
    fprintf(' Finding ergodic distribution...\n');
    distribution = double(full(ergodicdist(sparse(statetrans),1,ergodic_tol)));
