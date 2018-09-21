function state_dist = find_stationary(dieprob,ytrans,yPdist,nyP,nyF,nx,betatrans,...
                                N,xgrid_wide,sav_wide)
                        
% DISTRIBUTION

% adjust income transition probabilities by probability of death
ytrans_deathadjust = (1-dieprob)*ytrans + dieprob*repmat(yPdist',[nyP*nyF 1]);
% transition probabilities associated with (beta, yP, yF)
Pi_beta_yP_yF = kron(betatrans,ytrans_deathadjust);
% transition probabilities assigned to exogenous states
exog_trans = kron(Pi_beta_yP_yF,ones(nx,1)); 
% xprime
% xp = (1+r)*repmat(sav,[1 nyT]) + netymat;

grid_probabilities = zeros(N,N/nx);
% Perform interpolation separately for each (yP,yF,nb) combination
for col = 1:N/nx
    fspace = fundef({'spli',xgrid_wide(:,col),0,1});
    grid_probabilities(:,col) = funbas(fspace,sav_wide(:,col));
end
Q  = dprod(exog_trans,grid_probabilities);

state_dist      = full(ergodicdist(sparse(Q)));
statedist_orig_SS = state_dist;
% statedist_SS      = permute(reshape(statedist_SS,[nx,nyP*nyF*nb]),[2 1]);
% statedist_SS  = kron(yTdist,statedist_SS);


end