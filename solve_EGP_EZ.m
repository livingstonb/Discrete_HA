function [AYdiff,model] = solve_EGP_EZ(beta,p,xgrid,sgrid,agrid_short,prefs,income,Iterating)

    agrid = repmat(agrid_short,p.nyP*p.nyF*p.nb,1);
    
    %% CONSTRUCT EXPECTATIONS MATRIX                                     
    betagrid = beta + prefs.betagrid0;
    
    if p.IterateBeta == 1 && p.Display == 1
        msg = sprintf(' %3.3f',betagrid);
        disp([' Trying betagrid =' msg])
    end
    
    % initial guess for consumption function, stacked state combinations
    % column vector of length p.nx * p.nyP * p.nyF * p.nb
    con = (1/p.R) * repmat(xgrid.full(:),p.nb,1) + extracon;
    
    % initial guess for value function
    V = con;
    
    % discount factor matrix, 
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    betastacked = kron(betagrid,ones(p.nyP*p.nyF*p.nx,1));
    betastacked = sparse(diag(betastacked));

    % Expectations operator (conditional on yT)
    % square matrix of dim p.nx*p.nyP*p.nyF*p.nb
    Emat = kron(prefs.betatrans,kron(income.ytrans,speye(p.nx)));
    
    %% EGP Iteration
    iter = 1;
    cdiff = 1;
    while iter<p.max_iter && cdiff>p.tol_iter
        if iter==1
            conlast = con;
            Vlast   = V;
        else
            conlast = conupdate;
            Vlast   = Vupdate;
        end
        iter = iter + 1;
        
        % interpolate to get c(x') using c(x)
        
        % c(x) and V(x)
        conlast = reshape(conlast,[p.ns p.nyP p.nyF p.nb]);
        Vlast   = reshape(Vlast,[p.ns p.nyP p.nyF p.nb]);
        % c(x') and V(x')
        c_xp = zeros(p.ns,p.nyP,p.nyF,p.nb,p.nyT);
        V_xp = zeros(p.ns,p.nyP,p.nyF,p.nb,p.nyT);
        
        % x'(s)
        temp_sav = repmat(sgrid.full(:),p.nb,p.nyT);
        temp_inc = repmat(kron(income.netymat,ones(p.ns,1)),p.nb,1);
        xp_s = (1+p.r)*temp_sav + temp_inc;
        xp_s = reshape(xp_s,[p.ns p.nyP p.nyF p.nb p.nyT]);
        
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);
            coninterp = griddedInterpolant(xgrid.full(:,iyP,iyF),conlast(:,iyP,iyF,ib),'linear');
            c_xp(:,iyP,iyF,ib,:) = reshape(coninterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
            Vinterp = griddedInterpolant(xgrid.full(:,iyP,iyF),Vlast(:,iyP,iyF,ib),'linear');
            V_xp(:,iyP,iyF,ib,:) = reshape(Vinterp(xp_s_ib_iyF_iyP(:)),[],1,1,1,p.nyT);
        end
        end
        end
        
         % reshape to take expecation over yT first
        c_xp        = reshape(c_xp,[],p.nyT);
        xp_s	    = reshape(xp_s,[],p.nyT);

        % matrix of next period muc, muc(x',yP',yF)
        mucnext     = 
        
        
    end
    
end