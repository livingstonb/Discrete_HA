function modelupdate = find_stationary_adist(p,model,income,heterogeneity,grids)
    % Finds the stationary distribution and transition matrix for a given
    % grids.a.vec
    
    modelupdate = model;

    fprintf(' Computing state-to-state transition probabilities... \n');

    nx = size(grids.a.vec,1);
    if nx == p.nx
        netymat = income.netymatHJB;
    elseif nx == p.nx_KFE
        netymat = income.netymatKFE;
    end

    if numel(p.r) > 1
        r_col = kron(p.r',ones(nx*p.nyP*p.nyF,1));
        r_mat = reshape(r_col,[nx,p.nyP,p.nyF,numel(p.r)]);
    else
        r_mat = p.r;
    end

    % cash-on-hand as function of (a,yP,yF,yT)
    x = squeeze(grids.a.matrix(:,:,:,1)) + netymat;
    
    % saving interpolated onto this grid
    sav = zeros(nx,p.nyP,p.nyF,p.nb,p.nyT);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        x_iyP_iyF = x(:,iyP,iyF,:);
        sav_iyP_iyF_ib = model.savinterp{iyP,iyF,ib}(x_iyP_iyF(:));
        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_ib,[nx 1 1 1 p.nyT]);
    end
    end
    end
    
    aprime_live = (1+repmat(r_mat,[1,1,1,1,p.nyT])) .* sav;

    % transition matrix over (x,yP,yF,beta) full asset space
    modelupdate.statetrans = sparse(nx*p.nyP*p.nyF*p.nb,nx*p.nyP*p.nyF*p.nb);
    % create spline object
    fspace = fundef({'spli',grids.a.vec,0,1});
    % get interpolated probabilities and take expectation over yT
    interp_live = funbas(fspace,aprime_live(:));
    interp_live = reshape(interp_live,nx*p.nyP*p.nyF*p.nb,nx*p.nyT);
    interp_live = interp_live * kron(speye(nx),income.yTdist);
    if p.Bequests == 1
        interp_death = interp_live;
    else
        interp_death = sparse(nx*p.nyP*p.nyF*p.nb,nx);
        interp_death(:,1) = 1;
    end

    col = 1;
    for ib2 = 1:p.nb
    for iyF2 = 1:p.nyF
    for iyP2 = 1:p.nyP
        transcol_live = kron(income.ytrans_live(:,col),ones(nx,1));
        transcol_death = kron(income.ytrans_death(:,col),ones(nx,1));
        
        transcol_live = transcol_live .* interp_live;
        transcol_death = transcol_death .* interp_death;

        % add new column to transition matrix
        modelupdate.statetrans(:,nx*(col-1)+1:nx*col) = ...
            (1-p.dieprob)*transcol_live + p.dieprob*transcol_death;
        col = col + 1;
    end
    end
    end
    clear transcol_live transcol_death interp_live interp_death

    % stationary distribution over states
    fprintf(' Finding ergodic distribution...\n');
    
    q=zeros(1,nx*p.nyP*p.nyF*p.nb);
    % Create valid initial distribution for both yF & beta
    % Repmat automatically puts equal weight on each beta
    q(1,1:nx*p.nyP:end)=repmat(income.yFdist,p.nb,1) / p.nb;
    diff=1; 
    iter = 1;
    while diff>1e-8 && iter < 5e4
        z = q*modelupdate.statetrans;
        diff = norm(z-q);
        q = z;
        
        fprintf('  Diff = %5.3E, Iteration = %u \n',diff,iter);
        iter = iter + 1;
    end
    if iter >= 5e4
        error('No conv to statdist, diff = %5.3e',diff)
    end
    
%     [q,~] = eigs(modelupdate.statetrans',[],1,1);
%     q = q / sum(q(:));

    modelupdate.adiff = diff;
%     modelupdate.adist = reshape(full(q),[nx,p.nyP,p.nyF,p.nb]);
%     modelupdate.adiff = 1e-8;
    modelupdate.adist = reshape(full(q'),[nx,p.nyP,p.nyF,p.nb]);
    
    % get distribution over (x,yP,yF,beta)
    xdist = kron(income.yTdist,reshape(modelupdate.adist,nx,[]));
    modelupdate.xdist = reshape(xdist,[nx*p.nyT p.nyP p.nyF p.nb]);
    
    % Extend xvals to (nx*p.nyT,p.nyP,p.nyF,p.nyT)
    incvals = reshape(income.ymat,[p.nyP*p.nyF p.nyT]);
    incvals = permute(incvals,[2 1]);
    incvals = kron(incvals,ones(nx,1));
    incvals = reshape(incvals,[nx*p.nyT p.nyP p.nyF]);
    modelupdate.y_x = repmat(incvals,[1 1 1 p.nb]);
    modelupdate.nety_x = income.lumptransfer + (1-p.labtaxlow)*incvals - p.labtaxhigh*max(incvals-income.labtaxthresh,0);
    modelupdate.nety_x = repmat(modelupdate.nety_x,[1 1 1 p.nb]);
    modelupdate.xvals = repmat(grids.a.vec,[p.nyT p.nyP p.nyF p.nb]) + modelupdate.nety_x;
end
