function state_dist = find_stationary(dieprob,ytrans,yPdist,nyP,nyF,nx,betatrans,...
                                N,xgrid_wide,sav_wide,sgrid,borrow_lim,...
                                xmax,r,nyT,netymat,yTdist,yPtrans)
                        
% DISTRIBUTION

% transition probabilities associated with (beta, yP, yF)
Pi_beta_yP_yF = kron(betatrans,ytrans); 

% nx by N/nx matrix for net income, conditional on yT
netymat_wideT = cell(nyT,1);
for iyT = 1:nyT
    netymat_wideT{iyT} = reshape(netymat(:,iyT),[nx N/nx]);
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

Q = grid_probabilities;
state_dist      = full(ergodicdist(sparse(Q)));
statedist_orig_SS = state_dist;

end