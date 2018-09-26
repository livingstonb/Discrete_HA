function [xgridm,savm,avg_mpc] = mpc_forward(xgrid,p,income,sav,...
    betatrans)
ergodic_tol = 1e-7;
netymatm = reshape(income.netymat,[1 p.nyP p.nyF p.nyT]);

% create new, evenly spaced grid
nn = p.nxmpc;
NN = nn * p.nyP * p.nyF * p.nb;
xgridm = linspace(0,1,p.nxmpc)';
xgridm = p.borrow_lim + (p.xmax-p.borrow_lim)*xgridm;
gridspace = xgridm(2) - xgridm(1);

% get interpolated saving function on grid
savm = zeros(nn,p.nyP,p.nyF,p.nb);
for ib = 1:p.nb
for iyF = 1:p.nyF
for iyP = 1:p.nyP
    savinterp{iyP,iyF,ib} = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),sav.orig_wide(:,iyP,iyF,ib),'linear','linear');
    savm(:,iyP,iyF,ib) = savinterp{iyP,iyF,ib}(xgridm);
end
end
end

% create transition probability matrix
trans = kron(betatrans,kron(eye(p.nyF),income.yPtrans));
fspace = fundef({'spli',xgridm,0,1});
grid_probabilities = zeros(NN,NN);
col = 1;
for ib2 = 1:p.nb
for iyF2 = 1:p.nyF
for iyP2 = 1:p.nyP
    % xprime if no death
    xp = (1+p.r)*repmat(savm(:),p.nyT,1) + ...
        kron(squeeze(netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
    % xprime if death (i.e. saving in past period set to 0)
    xp_death = kron(squeeze(netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
    
    interp = (1-p.dieprob) * funbas(fspace,xp) + p.dieprob * funbas(fspace,xp_death);
    interp = reshape(full(interp),[],p.nyT*nn);
    % Multiply by yT distribution
    newcolumn = interp * kron(speye(nn),income.yTdist);
    % Multiply by (beta,yF,yP) distribution
    newcolumn = bsxfun(@times,kron(trans(:,col),ones(nn,1)),newcolumn);
    
    grid_probabilities(:,nn*(col-1)+1:nn*col) = newcolumn;
    col = col + 1;
end
end
end

% find stationary distribution
state_dist = full(ergodicdist(sparse(grid_probabilities),1,ergodic_tol));

% shift stationary distribution to get new grid
nn = nn + 1;
NN = nn * p.nyP * p.nyF * p.nb;
state_dist_wide = reshape(state_dist,[nn-1 p.nyP p.nyF p.nb]);
state_dist = zeros(nn,p.nyP,p.nyF,p.nb);
for ib = 1:p.nb
for iyF = 1:p.nyF
for iyP = 1:p.nyP
    state_dist(2:end,iyP,iyF,ib) = state_dist_wide(:,iyP,iyF,ib);
    state_dist(1,iyP,iyF,ib) = 0;
end
end
end
state_dist = state_dist(:);
xgridm(end+1) = xgridm(end) + gridspace;

% get interpolated saving function on new grid
savm = zeros(nn,p.nyP,p.nyF,p.nb);
for ib = 1:p.nb
for iyF = 1:p.nyF
for iyP = 1:p.nyP
    savinterp{iyP,iyF,ib} = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),sav.orig_wide(:,iyP,iyF,ib),'linear','linear');
    savm(:,iyP,iyF,ib) = savinterp{iyP,iyF,ib}(xgridm);
end
end
end

% find new transition probability matrix
fspace = fundef({'spli',xgridm,0,1});
grid_probabilities = zeros(NN,NN);
col = 1;
for ib2 = 1:p.nb
for iyF2 = 1:p.nyF
for iyP2 = 1:p.nyP
    % xprime if no death
    xp = (1+p.r)*repmat(savm(:),p.nyT,1) + ...
        kron(squeeze(netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
    % xprime if death (i.e. saving in past period set to 0)
    xp_death = kron(squeeze(netymatm(1,iyP2,iyF2,:)),ones(nn*p.nyP*p.nyF*p.nb,1));
    
    interp = (1-p.dieprob) * funbas(fspace,xp) + p.dieprob * funbas(fspace,xp_death);
    interp = reshape(full(interp),[],p.nyT*nn);
    % Multiply by yT distribution
    newcolumn = interp * kron(speye(nn),income.yTdist);
    % Multiply by (beta,yF,yP) distribution
    newcolumn = bsxfun(@times,kron(trans(:,col),ones(nn,1)),newcolumn);
    
    grid_probabilities(:,nn*(col-1)+1:nn*col) = newcolumn;
    col = col + 1;
end
end
end


% find new stationary distribution
SS_dist = full(ergodicdist(sparse(grid_probabilities),1,ergodic_tol));

% find future average consumption
savm = zeros(nn,p.nyP,p.nyF,p.nb);
for ib = 1:p.nb
for iyF = 1:p.nyF
for iyP = 1:p.nyP
    savm(:,iyP,iyF,ib) = savinterp{iyP,iyF,ib}(xgridm);
end
end
end
savm(savm<p.borrow_lim) = p.borrow_lim;
conm = repmat(xgridm,[1 p.nyP p.nyF p.nb]) - savm - p.savtax * max(savm - p.savtaxthresh,0);

% find mpcs

dist{1} = state_dist';
dist{2} = dist{1} * grid_probabilities;
dist{3} = dist{2} * grid_probabilities;
dist{4} = dist{3} * grid_probabilities;

% Later, use multiples of gridspace%
mpcamount = gridspace;
avg_mpc{1} = (dist{1} - SS_dist') * conm(:)/mpcamount;

for period = 2:4
    avg_mpc{period} = avg_mpc{period-1} + (dist{period} - SS_dist') * conm(:)/mpcamount;
end


end