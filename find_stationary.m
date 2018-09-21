function state_dist = find_stationary(dieprob,ytrans,yPdist,nyP,nyF,nx,betatrans,...
                                N,xgrid_wide,sav_wide,sgrid,borrow_lim,...
                                xmax,r,nyT,netymat,yTdist)
                        
% DISTRIBUTION

% adjust income transition probabilities by probability of death
ytrans_deathadjust = (1-dieprob)*ytrans + dieprob*repmat(yPdist',[nyP*nyF 1]);
% transition probabilities associated with (beta, yP, yF)
Pi_beta_yP_yF = kron(betatrans,ytrans_deathadjust);
% transition probabilities assigned to exogenous states
exog_trans = kron(Pi_beta_yP_yF,ones(nx,1)); 

% nx by N/nx matrix for net income, conditional on yT
netymat_wideT = cell(nyT,1);
for iyT = 1:nyT
    netymat_wideT{iyT} = reshape(netymat(:,iyT),[nx N/nx]);
end

for col2 = 1:N/nx
    fspace = fundef({'spli',xgrid_wide(:,col2),0,1});
    for col1 = 1:N/nx
        col1_col2_probs = 0;
        for iyT = 1:nyT
             xp = (1+r)*repmat(sav_wide(:,col1),[1 nyT]) + netymat_wideT{iyT}(:,col2);
             col1_col2_probs = col1_col2_probs + yTdist(iyT) * funbas(fspace,xp) .* Pi_beta_yP_yF(col1,col2) ;
        end
        grid_probabilities(nx*(col1-1)+1:nx*col1,nx*(col2-1)+1:nx*col2) = col1_col2_probs;
    end
end

Q = grid_probabilities;
state_dist      = full(ergodicdist(sparse(Q)));
statedist_orig_SS = state_dist;
% statedist_SS      = permute(reshape(statedist_SS,[nx,nyP*nyF*nb]),[2 1]);
% statedist_SS  = kron(yTdist,statedist_SS);


end