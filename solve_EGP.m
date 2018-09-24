function [AYdiff,con_opt,sav_opt,state_dist,cdiff] = solve_EGP(beta,p,...
    xgrid_wide,ytrans,betatrans,sgrid_wide,u1,u1inv,netymat,...
    yTdist,beq1,yPtrans,meany,ergodic_tol)

nx = p.nx;
ns = p.ns;
nyF = p.nyF;
nyP = p.nyP;
nyT = p.nyT;
nb = p.nb;
N = p.N;
r = p.r;
targetAY = p.targetAY;
max_iter = p.max_iter;
tol_iter = p.tol_iter;
dieprob = p.dieprob;
savtax = p.savtax;
savtaxthresh = p.savtaxthresh;
borrow_lim = p.borrow_lim;
temptation = p.temptation;
betawidth = p.betawidth;


if  nb == 1
    betagrid = beta;
elseif nb ==2 
    betagrid = [beta-betawidth;beta+betawidth];
end

if p.IterateBeta == 1
    disp(['Trying betagrid = ' num2str(betagrid)])
end

% initial guess for consumption function
con = r * xgrid_wide(:);

% discount factor matrix
betastacked = reshape(repmat(betagrid',nyP*nyF*nx,1),N,1);

% Expectations operator (conditional on yT)
Emat = kron(betatrans,kron(ytrans,speye(nx)));

iter = 1;
cdiff = 1;
%% EGP ITERATION
while iter<max_iter && cdiff>tol_iter
    if iter==1
        conlast = con;
    else
        conlast = conupdate;
    end
    iter = iter + 1;

    % interpolate to get c(x') using c(x)
    conlast_wide = reshape(conlast,ns,nyP*nyF*nb);
    % initialize cons as function of x',yT
    c_xp = zeros(N,nyT);
    
    x_s = (1+r)*repmat(sgrid_wide(:),1,nyT) + netymat;
    
    % c(x')
    c_xp = zeros(ns,nyP*nyF*nb,nyT);
    for col = 1:nyP*nyF*nb
        coninterp = griddedInterpolant(xgrid_wide(:,col),conlast_wide(:,col),'linear','linear');
        for iyT = 1:nyT
            x_s_wide = reshape(x_s(:,iyT),ns,nyP*nyF*nb); 
            c_xp(:,col,iyT) = coninterp(x_s_wide(:,col));
        end
    end
    c_xp = reshape(c_xp,[],nyT);
    
    mucnext  = u1(c_xp) - temptation/(1+temptation)*u1(x_s);
    % muc this period as a function of s
    muc_s = (1-dieprob)*(1+r)*betastacked.*(Emat*(mucnext*yTdist))./(1+savtax.*(sgrid_wide(:)>=savtaxthresh))...
        + dieprob*beq1(sgrid_wide(:));
    % _wide variables have dimension nx by nyP*nyF*nb, or nx by N/nx

    % consumption as a function of s
    con_s = u1inv(muc_s);
    % cash-in-hand (x) as a function of s
    x_s = sgrid_wide(:) + savtax * max(sgrid_wide(:)-savtaxthresh,0) + con_s;

    % interpolate from x(s) to get s(x), interpolate for each (beta,yP,yF)
    % separately
    x_s_wide = reshape(x_s,ns,N/ns);
    sav_wide = zeros(ns,N/ns);
    for col=1:N/ns
        sav_wide(:,col) = interp1(x_s_wide(:,col),sgrid_wide(:,col),xgrid_wide(:,col),'linear','extrap'); 
    end


    % deal with borrowing limit
    sav_wide(sav_wide<borrow_lim) = borrow_lim;
    sav_opt = sav_wide(:);

    conupdate = xgrid_wide(:) - sav_opt - savtax * max(sav_opt-savtaxthresh,0);

    cdiff = max(abs(conupdate-conlast));
    if mod(iter,50) ==0
    	disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
    end
end

con_opt = conupdate;

%% DISTRIBUTION
fprintf(' Computing state-to-state transition probabilities... \n');

% Use original xgrids
yFtrans = eye(nyF);
xgridm = reshape(xgrid_wide,[nx nyP nyF nb]);
savm = reshape(sav_opt,[nx nyP nyF nb]);
netymatm = reshape(netymat,[nx nyP nyF nb nyT]);

grid_probabilities = zeros(N,N);
outerblock = 1;
for ib2 = 1:nb
for iyF2 = 1:nyF
for iyP2 = 1:nyP
    fspace = fundef({'spli',xgridm(:,iyP2,iyF2,ib2),0,1});
    
    newcolumn = zeros(N,nx);
    innerblock = 1;
    for ib1 = 1:nb
    for iyF1 = 1:nyF
    for iyP1 = 1:nyP
        state1prob = 0;
        for iyT = 1:nyT
            xp = (1+r)*savm(:,iyP1,iyF1,ib1) + netymatm(:,iyP2,iyF2,ib2,iyT);
            state1prob = state1prob + yTdist(iyT) * yPtrans(iyP1,iyP2) * yFtrans(iyF1,iyF2) * betatrans(ib1,ib2) * funbas(fspace,xp);
        end
        newcolumn(nx*(innerblock-1)+1:nx*innerblock,:) = state1prob;
        innerblock = innerblock + 1;
    end
    end
    end
    grid_probabilities(:,nx*(outerblock-1)+1:nx*outerblock) = newcolumn;
    outerblock = outerblock + 1;
end
end
end

% SS probability of residing in each state
fprintf(' Finding ergodic distribution...\n');
state_dist      = full(ergodicdist(sparse(grid_probabilities),1,ergodic_tol));

% SS wealth/gross income ratio
mean_s = sav_opt' * state_dist;
fprintf(' A/Y = %2.3f\n',mean_s/meany);
%AYdiffsq = (mean_s/meany - targetAY)^2;
AYdiff = mean_s/meany -  targetAY;

end