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
% xprime
xp = (1+r)*repmat(sav_wide(:),[1 nyT]) + netymat;

% Perform interpolation separately for each (yP,yF,nb,yT) combination
grid_probabilities = zeros(N,N);

% Initial (yP,yF,beta)
for col = 1:N/nx
fspace = fundef({'spli',xgrid_wide(:,col),0,1});
    % Next (yP,yF,beta)
    for col2 = 1:N/nx
        for iyT = 1:nyT
            xp_wide = reshape(xp(:,iyT),nx,N/nx);
            grid_probabilities(nx*(col-1)+1:nx*col,nx*(col2-1)+1:nx*col2) = grid_probabilities(nx*(col-1)+1:nx*col,nx*(col2-1)+1:nx*col2) + yTdist(iyT) * funbas(fspace,xp_wide(:,col)) .* exog_trans(nx*(col-1)+1:nx*col,col2) ;
        end
    end
end

Q = grid_probabilities;
state_dist      = full(ergodicdist(sparse(Q)));
statedist_orig_SS = state_dist;
% statedist_SS      = permute(reshape(statedist_SS,[nx,nyP*nyF*nb]),[2 1]);
% statedist_SS  = kron(yTdist,statedist_SS);


end