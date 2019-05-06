function income = gen_income_variables(p,prefs)
    % Given a structure of parameters, p, this function generates a
    % structure variable called 'income' with fields associated with the
    % specified income distribution
    
    LoadIncome = ~isempty(p.IncomeProcess); % 1 if load from file
    
    if LoadIncome==1
        Import = load(p.IncomeProcess);
    end

    %% PERSISTENT INCOME
    % rowenhurst
    if LoadIncome == 1
        logyPgrid = Import.discmodel1.logyPgrid;
        yPdist = Import.discmodel1.yPdist;
        yPtrans = Import.discmodel1.yPtrans;
        p.nyP = length(logyPgrid);
        logyPgrid = reshape(logyPgrid,[],1);
        yPdist = reshape(yPdist,[],1);
    elseif p.nyP>1
        [logyPgrid, yPtrans, yPdist] = rouwenhorst(p.nyP, -0.5*p.sd_logyP^2, p.sd_logyP, p.rho_logyP);
    else
        logyPgrid = 0;
        yPdist = 1;
        yPtrans = 1;
    end  

    yPgrid = exp(logyPgrid);
    yPgrid = yPgrid/(yPdist'*yPgrid);
    logyPgrid = log(yPgrid);
    yPcumdist = cumsum(yPdist,1);
    yPcumtrans = cumsum(yPtrans,2);
    
    %% TRANSITORY INCOME
    % disretize normal distribution
    if LoadIncome == 1
        logyTgrid = Import.discmodel1.logyTgrid;
        yTdist = Import.discmodel1.yTdist;
        p.nyT = length(logyTgrid);
        logyTgrid = reshape(logyTgrid,[],1);
        yTdist = reshape(yTdist,[],1);
        yTcumdist = cumsum(yTdist,1);
    elseif p.nyT>1

        %moments of mixture distribution
        lmu2 = p.lambdaT.*p.sd_logyT^2;
        lmu4 = 3.*p.lambdaT.*(p.sd_logyT^4);

        %fit thjose moments
        optionsNLLS = optimoptions(@lsqnonlin,'Display','Off');
        lpar = lsqnonlin(@(lp)discretize_normal_var_kurt(lp,p.nyT,-lmu2/2,lmu2,lmu4),[2 0.1],[],[],optionsNLLS);
        [lf,lx,lp] = discretize_normal_var_kurt(lpar,p.nyT,-lmu2/2,lmu2,lmu4);
        logyTgrid = lx;
        yTdist = lp;
        yTcumdist = cumsum(yTdist,1);

    elseif p.nyT==1
        logyTgrid = 0;
        yTdist = 1;
        yTcumdist = 1;
    end
    
    yTgrid = exp(logyTgrid);
    yTgrid = yTgrid/(yTdist'*yTgrid);
    logyTgrid = log(yTgrid);

    %% FIXED EFFECT
    if p.nyF>1
        width = fzero(@(x)discrete_normal(p.nyF,-0.5*p.sd_logyF^2 ,p.sd_logyF ,x),2);
        [~,logyFgrid,yFdist] = discrete_normal(p.nyF,-0.5*p.sd_logyF^2 ,p.sd_logyF ,width);
        logyFgrid = reshape(logyFgrid,[],1);
        yFdist = reshape(yFdist,[],1);
    elseif p.nyF==1
        logyFgrid = 0;
        yFdist = 1;
    end
    yFgrid = exp(logyFgrid);
    % normalize fixed effect such that mean = 1 if annual, 1/4 if quarterly
    yFgrid = yFgrid/(yFdist'*yFgrid*p.freq);
    logyFgrid = log(yFgrid);
    yFcumdist = cumsum(yFdist,1);
    
    if size(yTgrid,2)>1 || size(yFgrid,2)>1 || size(yPgrid,2)>1
        error('All income grids must be column vectors')
    end
    if size(yTdist,2)>1 || size(yFdist,2)>1 || size(yPdist,2)>1
        error('All income distributions must be column vectors')
    end

    %% OTHER INCOME VARIABLES
    % transition probabilities for yP-yF combined grid
    ytrans = kron(eye(p.nyF),yPtrans);

    % construct matrix of y combinations
    ymat = repmat(yPgrid,p.nyF,1) .* kron(yFgrid,ones(p.nyP,1)) * yTgrid';

    % distribution of ymat
    ymatdist = repmat(yPdist,p.nyF,1) .* kron(yFdist,ones(p.nyP,1)) * yTdist';

    % find mean y
    % isolate unique (yT,yF,yP) combinations
    temp = sortrows([ymat(:) ymatdist(:)],1);
    ysort = temp(:,1);
    ysortdist = temp(:,2);
    ycumdist_sort = cumsum(ysortdist);
    
    % 1-period statistics
    meany1 = ymat(:)'*ymatdist(:);
    totgrossy1 = meany1;

    % find tax threshold on labor income
    if numel(ysort)>1
        labtaxthresh = lininterp1(ycumdist_sort,ysort,p.labtaxthreshpc);
    else
        labtaxthresh = 0;
    end    

    % find net income
    totgrossyhigh = max(ymat(:)-labtaxthresh,0)'*ymatdist(:);
    lumptransfer = p.labtaxlow*totgrossy1 + p.labtaxhigh*totgrossyhigh;
    % netymat is N by nyT matrix
    netymat = lumptransfer + (1-p.labtaxlow)*ymat - p.labtaxhigh*max(ymat-labtaxthresh,0);
    meannety1 = netymat(:)'*ymatdist(:);
    
    % net y values on HJB grid
    netymat_temp = reshape(netymat,[1 p.nyP p.nyF p.nyT]);
    netymatHJB = repmat(netymat_temp,[p.nx 1 1 1]);
    netymatKFE = repmat(netymat_temp,[p.nx_KFE 1 1 1]);

    % full transition matrix with beta and IES transitions, excluding and including death
    if p.ResetIncomeUponDeath == 1
        yPtrans_death = repmat(yPdist',p.nyP,1);
    else
        yPtrans_death = yPtrans;
    end

    if (numel(p.risk_aver) == 1) && (numel(p.invies) == 1)
        ytrans_live = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans));
        ytrans_death = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans_death));
    else
        ytrans_live = kron(prefs.IEStrans,kron(eye(p.nyF),yPtrans));
        ytrans_death = kron(prefs.IEStrans,kron(eye(p.nyF),yPtrans_death));
    end
    
        % Store income variables in a structure
    newfields = {'ymat','netymat','meany1','yPgrid',...
        'yTgrid','yFgrid','yPdist','yTdist','yFdist','yPcumtrans',...
        'yPtrans','yPcumdist','yFcumdist','yTcumdist','ytrans',...
        'meannety1','labtaxthresh','lumptransfer','ysortdist','ysort',...
        'ymatdist','netymatHJB','netymatKFE','ytrans_live','ytrans_death'};
    for i = 1:numel(newfields)
        income.(newfields{i}) = eval(newfields{i});
    end



end