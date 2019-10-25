function [AYdiff,modelupdate] = find_stationary_adist(p,model,income,grids,heterogeneity)
    % Finds the stationary distribution and transition matrix for a given
    % grids.a.vec
    
    %% ----------------------------------------------------------------
    % FIND STATIONARY DISTRIBUTION
    % -----------------------------------------------------------------
    modelupdate = model;

    fprintf(' Computing state-to-state transition probabilities... \n');

    nx = size(grids.a.vec,1);
    if nx == p.nx
        netymat = income.netymatEGP;
    elseif nx == p.nx_DST
        netymat = income.netymatDST;
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
        x_iyP_iyF_ib = x(:,iyP,iyF,ib,:);
        sav_iyP_iyF_ib = model.savinterp{iyP,iyF,ib}(x_iyP_iyF_ib(:));
        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_ib,[nx 1 1 1 p.nyT]);
    end
    end
    end

    % transition matrix over (x,yP,yF,beta) full asset space
    modelupdate.statetrans = get_transition_matrix(p,income,grids,nx,sav,r_mat);

    % stationary distribution over states
    fprintf(' Finding ergodic distribution...\n');
    q = get_distribution(p,income,nx,modelupdate.statetrans,heterogeneity);
%     [q,~] = eigs(modelupdate.statetrans',[],1,1);
%     q = q / sum(q(:));

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

    %% ----------------------------------------------------------------
    % POLICY FUNCTIONS ETC...
    % -----------------------------------------------------------------
    % get saving policy function defined on xgrid
    modelupdate.sav_x = zeros(p.nx_DST*p.nyT,p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP 
        modelupdate.sav_x(:,iyP,iyF,ib) = modelupdate.savinterp{iyP,iyF,ib}(modelupdate.xvals(:,iyP,iyF,ib));
    end
    end
    end
    modelupdate.sav_x = max(modelupdate.sav_x,p.borrow_lim);

    % Collapse the asset distribution from (a,yP_lag,yF_lag,beta_lag) to (a,beta_lag) for norisk
    % model, and from (x,yP,yF,beta) to (x,beta)
    if p.nyP>1 && p.nyF>1
        % a
        modelupdate.adist_noincrisk =  sum(sum(modelupdate.adist,3),2);
        % x
        modelupdate.xdist_noincrisk    = sum(sum(modelupdate.xdist,3),2);
    elseif (p.nyP>1 && p.nyF==1) || (p.nyP==1 && p.nyF>1)
        modelupdate.adist_noincrisk =  sum(modelupdate.adist,2);
        modelupdate.xdist_noincrisk    = sum(modelupdate.xdist,2);
    elseif p.nyP==1 && p.nyF==1
        modelupdate.adist_noincrisk = modelupdate.adist;
        modelupdate.xdist_noincrisk    = modelupdate.xdist;
    end

    % Policy functions associated with xdist
    modelupdate.con_x = modelupdate.xvals - modelupdate.sav_x ...
    	- p.savtax*max(modelupdate.sav_x-p.savtaxthresh,0);
    
    % mean saving, mean assets
	modelupdate.mean_a = modelupdate.adist(:)' * grids.a.matrix(:);
    
    mean_assets = modelupdate.mean_a;
    fprintf(' A/Y = %2.5f\n',mean_assets/(income.meany1*p.freq));
    AYdiff = mean_assets/(income.meany1*p.freq) -  p.targetAY;
end

%% ----------------------------------------------------------------
% TRANSITION MATRIX
% -----------------------------------------------------------------
function trans = get_transition_matrix(p,income,grids,nx,sav,r_mat)
	aprime_live = (1+repmat(r_mat,[1,1,1,1,p.nyT])) .* sav;

	% create interpolant object
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

    trans = sparse(nx*p.nyP*p.nyF*p.nb,nx*p.nyP*p.nyF*p.nb);
    col = 1;
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        transcol_live = kron(income.ytrans_live(:,col),ones(nx,1));
        transcol_death = kron(income.ytrans_death(:,col),ones(nx,1));
        
        transcol_live = transcol_live .* interp_live;
        transcol_death = transcol_death .* interp_death;

        % add new column to transition matrix
        trans(:,nx*(col-1)+1:nx*col) = ...
            (1-p.dieprob)*transcol_live + p.dieprob*transcol_death;
        col = col + 1;
    end
    end
    end
end

%% ----------------------------------------------------------------
% iTERATIVE METHOD TO FIND STATIONARY DISTRIBUTION
% -----------------------------------------------------------------
function q = get_distribution(p,income,nx,statetrans,heterogeneity)
	q=ones(nx,p.nyP,p.nyF,p.nb);

    % create valid initial distribution for both yF & beta
    beta_dist = reshape(heterogeneity.betadist,[1,1,1,p.nbeta]);
    yFdist = reshape(income.yFdist,[1,1,p.nyF,1]);
    q = beta_dist .* yFdist;

    if (p.nbeta==1) && (p.nb>1)
        q = repmat(q,[nx,p.nyP,1,p.nb]);
    else
        q = repmat(q,[nx,p.nyP,1,1]);
    end
    q = q(:)' / sum(q(:));

    diff = 1; 
    iter = 1;
    while diff>1e-9 && iter < 5e4
        z = q * statetrans;
        diff = norm(z-q);
        q = z;
        
        if mod(iter,50) == 0
            fprintf('  Diff = %5.3E, Iteration = %u \n',diff,iter);
        end
        iter = iter + 1;
    end

    if iter >= 5e4
        error('No conv to statdist, diff = %5.3e',diff)
    end
end