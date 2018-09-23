function [AYdiffsq,con_opt,sav_opt,state_dist,cdiff] = solve_EGP(beta,p,...
    xgrid_wide,ytrans,betatrans,sgrid_wide,u1,u1inv,netymat,meany,...
    yTdist,beq1)

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


if  nb == 1
    betagrid = beta;
elseif nb ==2 
    betagrid = [beta-betawidth;beta+betawidth];
end

disp(['Trying betagrid = ' num2str(betagrid)])

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
    % get consumption as function of x' rather than s
    for iyT = 1:nyT
        x_s_wide = reshape(x_s(:,iyT),ns,nyP*nyF*nb); 
        c_xpT_wide = zeros(ns,nyP*nyF*nb);
        for col = 1:nyP*nyF*nb
            interpolant = griddedInterpolant(xgrid_wide(:,col),conlast_wide(:,col),'linear','nearest');
            c_xpT_wide(:,col) = interpolant(x_s_wide(:,col));
            %c_xpT_wide(:,col) = interp1(xgrid_wide(:,col),conlast_wide(:,col),x_s_wide(:,col),'linear','extrap');
        end
        c_xp(:,iyT)  = c_xpT_wide(:);
    end

    mucnext  = u1(c_xp);
    % muc this period as a function of s
    muc_s = (1-dieprob)*(1+r)*betastacked.*(Emat*(mucnext*yTdist))...
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
    constr = bsxfun(@lt,xgrid_wide,x_s_wide(1,:));
    sav_wide(constr) = borrow_lim;
    sav_opt = sav_wide(:);

    conupdate = xgrid_wide(:) - sav_opt - savtax * max(sav_opt-savtaxthresh,0);

    cdiff = max(abs(conupdate-conlast));
    if mod(iter,50) ==0
    	disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
    end
end

con_opt = conupdate;

%% DISTRIBUTION
% transition probabilities associated with (beta, yP, yF)
Pi_beta_yP_yF = kron(betatrans,ytrans); 

% nx by N/nx matrix for net income, conditional on yT
netymat_wideT = cell(nyT,1);
for iyT = 1:nyT
    netymat_wideT{iyT} = reshape(netymat(:,iyT),nx,N/nx);
end

for col2 = 1:N/nx
    fspace = fundef({'spli',xgrid_wide(:,col2),0,1});
    xmax_col2 = max(xgrid_wide(:,col2));
    xmin_col2 = min(xgrid_wide(:,col2));
    for col1 = 1:N/nx
        col1_col2_probs = 0;
        for iyT = 1:nyT
             xp = (1+r)*sav_wide(:,col1) + netymat_wideT{iyT}(:,col2);
             xp = min(max(xp,xmin_col2),xmax_col2);
             col1_col2_probs = col1_col2_probs + yTdist(iyT) * funbas(fspace,xp) .* Pi_beta_yP_yF(col1,col2);
        end
        grid_probabilities(nx*(col1-1)+1:nx*col1,nx*(col2-1)+1:nx*col2) = col1_col2_probs;
    end
end

% SS probability of residing in each state
state_dist      = full(ergodicdist(sparse(grid_probabilities)));

% SS wealth/gross income ratio
AYdiffsq = (sav_opt' * state_dist(:)/meany - targetAY)^2;

end