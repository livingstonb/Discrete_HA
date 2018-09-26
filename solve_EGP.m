function [AYdiff,con_opt,sav_opt,conm,savm,state_dist,cdiff,xgridm] = solve_EGP(beta,p,...
    xgrid,sgrid,betatrans,u1,beq1,u1inv,ergodic_tol,income,ExpandGrid)

ytrans = income.ytrans;
netymat = income.netymat;
yTdist = income.yTdist;
yPtrans = income.yPtrans;
meany = income.meany;

if  p.nb == 1
    betagrid = beta;
elseif nb ==2 
    betagrid = [beta-betawidth;beta+betawidth];
end

if p.IterateBeta == 1
    disp(['Trying betagrid = ' num2str(betagrid)])
end

% initial guess for consumption function
con = p.r * xgrid.orig_wide(:);

% discount factor matrix
betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));

% Expectations operator (conditional on yT)
Emat = kron(betatrans,kron(ytrans,speye(p.nx)));

iter = 1;
cdiff = 1;
%% EGP ITERATION
while iter<p.max_iter && cdiff>p.tol_iter
    if iter==1
        conlast = con;
    else
        conlast = conupdate;
    end
    iter = iter + 1;

    % interpolate to get c(x') using c(x)
    conlast_wide = reshape(conlast,[p.ns p.nyP p.nyF p.nb]);
    
    x_s_wide = (1+p.r)*repmat(sgrid.wide(:),1,p.nyT) + kron(netymat,ones(p.ns,1));
    x_s_wide = reshape(x_s_wide,[p.ns p.nyP p.nyF p.nb p.nyT]);
    
    % c(x')
    c_xp = zeros(p.ns,p.nyP,p.nyF,p.nb,p.nyT);
    for ib = 1:p.nb
    for iyF  = 1:p.nyF
    for iyP = 1:p.nyP
        coninterp = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),conlast_wide(:,iyP,iyF,ib),'linear','linear');
        c_xp_temp = coninterp(reshape(x_s_wide(:,iyP,iyF,ib,:),[],1));
        c_xp(:,iyP,iyF,ib,:) = reshape(c_xp_temp,[],1,1,1,p.nyT);
    end
    end
    end
    c_xp = reshape(c_xp,[],p.nyT);
    x_s_wide = reshape(x_s_wide,[],p.nyT);
    
    mucnext  = u1(c_xp) - p.temptation/(1+p.temptation)*u1(x_s_wide);
    % muc this period as a function of s
    muc_s = (1-p.dieprob)*(1+p.r)*betastacked.*(Emat*(mucnext*yTdist))./(1+p.savtax.*(sgrid.wide(:)>=p.savtaxthresh))...
        + p.dieprob*beq1(sgrid.wide(:));
    % _wide variables have dimension nx by nyP*nyF*nb, or nx by N/nx

    % consumption as a function of s
    con_s = u1inv(muc_s);
    % cash-in-hand (x) as a function of s
    x_s_wide = sgrid.wide(:) + p.savtax * max(sgrid.wide(:)-p.savtaxthresh,0) + con_s;

    % interpolate from x(s) to get s(x), interpolate for each (beta,yP,yF)
    % separately
    x_s_wide = reshape(x_s_wide,[p.ns p.nyP p.nyF p.nb]);
    sav_wide = zeros(p.ns,p.nyP,p.nyF,p.nb);
    sgrid_reshaped = reshape(sgrid.wide,[p.ns p.nyP p.nyF p.nb]);
    for ib = 1:p.nb
    for iyF  = 1:p.nyF
    for iyP = 1:p.nyP
        savinterp = griddedInterpolant(x_s_wide(:,iyP,iyF,ib),sgrid_reshaped(:,iyP,iyF,ib),'linear','linear');
        sav_wide(:,iyP,iyF,ib) = savinterp(xgrid.orig_wide(:,iyP,iyF,ib)); 
    end
    end
    end

    % deal with borrowing limit
    sav_wide(sav_wide<p.borrow_lim) = p.borrow_lim;
    sav_opt = sav_wide(:);

    conupdate = xgrid.orig_wide(:) - sav_opt - p.savtax * max(sav_opt-p.savtaxthresh,0);

    cdiff = max(abs(conupdate-conlast));
    if mod(iter,50) ==0
    	disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
    end
    
end

con_opt = conupdate;

%% DISTRIBUTION
fprintf(' Computing state-to-state transition probabilities... \n');
savm = reshape(sav_opt,[p.nx p.nyP p.nyF p.nb]);

% Create long grid
if ExpandGrid == 1
    xgridm = linspace(0,1,p.nxlong)';
    xgridm = repmat(xgridm,[1 p.nyP p.nyF]) .^(1/p.xgrid_par);
    netymatm = reshape(netymat,[1 p.nyP p.nyF p.nyT]);
    netymatm = repmat(netymatm,[p.nxlong 1 1 1]);
    xgridm = p.borrow_lim + min(netymatm,[],4) + (p.xmax-p.borrow_lim)*xgridm;

    savlong = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        savinterp = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),savm(:,iyP,iyF,ib),'linear','linear');
        savlong(:,iyP,iyF,ib) = savinterp(xgridm(:,iyP,iyF));
    end
    end
    end
    savm = savlong;
    NN = p.nxlong * p.nyP * p.nyF * p.nb;
    nn = p.nxlong;
else % Use original xgrids
    netymatm = reshape(netymat,[1 p.nyP p.nyF p.nyT]);
    netymatm = repmat(netymatm,[p.nx 1 1 1]);
    xgridm = xgrid.orig_wide;
    NN = p.N;
    nn = p.nx;
end

trans = kron(betatrans,kron(eye(p.nyF),yPtrans));
grid_probabilities = zeros(NN,NN);
col = 1;
for ib2 = 1:p.nb
for iyF2 = 1:p.nyF
for iyP2 = 1:p.nyP
    fspace = fundef({'spli',xgridm(:,iyP2,iyF2,ib2),0,1});
    % xprime if no death
    xp = (1+p.r)*repmat(savm(:),p.nyT,1) + ...
        kron(squeeze(netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
    % xprime if death (i.e. saving in past period set to 0)
    xp_death = kron(squeeze(netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
    
    interp = (1-p.dieprob) * funbas(fspace,xp) + p.dieprob * funbas(fspace,xp_death);
    interp = reshape(full(interp),[],p.nyT*nn);
    % Multiply by yT distribution
    newcolumn = interp * kron(speye(nn),yTdist);
    % Multiply by (beta,yF,yP) distribution
    newcolumn = bsxfun(@times,kron(trans(:,col),ones(nn,1)),newcolumn);
    
    grid_probabilities(:,nn*(col-1)+1:nn*col) = newcolumn;
    col = col + 1;
end
end
end

% SS probability of residing in each state
fprintf(' Finding ergodic distribution...\n');
state_dist      = full(ergodicdist(sparse(grid_probabilities),1,ergodic_tol));

% SS wealth/gross income ratio
mean_s = savm(:)' * state_dist;

% policy functions
savm = savm(:);
conm = xgridm(:) - savm(:) - p.savtax*max(savm(:)-p.savtaxthresh,0);

fprintf(' A/Y = %2.3f\n',mean_s/meany);
%AYdiffsq = (mean_s/meany - targetAY)^2;
AYdiff = mean_s/meany -  p.targetAY;

end