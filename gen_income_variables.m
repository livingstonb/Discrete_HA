function income = gen_income_variables(p)

    %% PERSISTENT INCOME
    % rowenhurst
    if p.LoadIncomeProcess == 1
        logyPgrid = load('QuarterlyIncomeDynamics/TransitoryContinuous/logyPgrid.txt');
        yPdist = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPdist.txt');
        yPtrans = load('QuarterlyIncomeDynamics/TransitoryContinuous/yPtrans.txt');
        p.nyP = length(logyPgrid);
        logyPgrid = logyPgrid';
    elseif p.nyP>1
        [logyPgrid, yPtrans, yPdist] = rouwenhorst(p.nyP, -0.5*p.sd_logyP^2, p.sd_logyP, p.rho_logyP);
    else
        logyPgrid = 0;
        yPdist = 1;
        yPtrans = 1;
    end  

    yPgrid = exp(logyPgrid);
    yPcumdist = cumsum(yPdist,1);
    yPcumtrans = cumsum(yPtrans,2);
    
    if size(yPgrid,2)>1
        error('yPgrid is a row vector, must be column vector')
    end
    if size(yPdist,2)>1
        error('yPdist is a row vector, must be column vector')
    end
    
    %% TRANSITORY INCOME
    % disretize normal distribution
    if p.LoadIncomeProcess == 1
        p.sig2T = load('QuarterlyIncomeDynamics/TransitoryContinuous/sig2T.txt');
        p.lambdaT = load('QuarterlyIncomeDynamics/TransitoryContinuous/lambdaT.txt');
    end

    if p.nyT>1

        %moments of mixture distribution
        lmu2 = p.lambdaT.*p.sd_logyT^2;
        lmu4 = 3.*p.lambdaT.*(p.sd_logyT^4);

        %fit thjose moments
        optionsNLLS = optimoptions(@lsqnonlin,'Display','Off');
        lpar = lsqnonlin(@(lp)discretize_normal_var_kurt(lp,p.nyT,lmu2,lmu4),[2 0.1],[],[],optionsNLLS);
        [lf,lx,lp] = discretize_normal_var_kurt(lpar,p.nyT,lmu2,lmu4);
        logyTgrid = lx;
        yTdist = lp;
        yTcumdist = cumsum(yTdist,1);

    elseif nyT==1
        logyTgrid = 0;
        yTdist = 1;
    end
    
    yTgrid = exp(logyTgrid);
    
    if size(yTgrid,2)>1
        error('yTgrid is a row vector, must be column vector')
    end
    if size(yTdist,2)>1
        error('yTdist is a row vector, must be column vector')
    end

    %% FIXED EFFECT
    if p.nyF>1
        width = fzero(@(x)discrete_normal(p.nyF,-0.5*p.sd_logyF^2 ,p.sd_logyF ,x),2);
        [~,logyFgrid,yFdist] = discrete_normal(p.nyF,-0.5*p.sd_logyF^2 ,p.sd_logyF ,width);
    elseif p.nyF==1
        logyFgrid = 0;
        yFdist = 1;
    end
    yFgrid = exp(logyFgrid);
    yFcumdist = cumsum(yFdist,1);

    if size(yFgrid,2)>1
        error('yFgrid is a row vector, must be column vector')
    end
    if size(yFdist,2)>1
        error('yFdist is a row vector, must be column vector')
    end

    %% OTHER INCOME VARIABLES
    % transition probabilities for yP-yF combined grid
    ytrans = kron(eye(p.nyF),yPtrans);

    % construct matrix of y combinationsx
    ymat = repmat(yPgrid,p.nyF,1) .* kron(yFgrid,ones(p.nyP,1)) * yTgrid';

    % distribution of ymat (values are repeated nb*nx times)
    ymatdist = repmat(yPdist,p.nyF,1) .* kron(yFdist,ones(p.nyP,1)) * yTdist';

    % find mean y
    % isolate unique (yT,yF,yP) combinations
    temp = sortrows([ymat(:) ymatdist(:)],1);
    ysort = temp(:,1);
    ysortdist = temp(:,2);
    ycumdist = cumsum(ysortdist);
    meany = ymat(:)'*ymatdist(:);
    original_meany = meany;

    % normalize so that mean annual gross income is 1
    if p.NormalizeY == 1
        ymat = ymat/(meany*p.freq);
        ysort = ysort/(meany*p.freq);
        meany = 1/p.freq;
    end
    totgrossy = meany;

    % find tax threshold on labor income
    if numel(ysort)>1
        labtaxthresh = lininterp1(ycumdist,ysort,p.labtaxthreshpc);
    else
        labtaxthresh = 0;
    end    

    % find net income
    totgrossyhigh = max(ymat(:)-labtaxthresh,0)'*ymatdist(:);
    lumptransfer = p.labtaxlow*totgrossy + p.labtaxhigh*totgrossyhigh;
    % netymat is N by nyT matrix
    netymat = lumptransfer + (1-p.labtaxlow)*ymat - p.labtaxhigh*max(ymat-labtaxthresh,0);
    meannety = netymat(:)'*ymatdist(:);
    
        % Store income variables in a structure
    newfields = {'ymat','netymat','meany','original_meany','yPgrid',...
        'yTgrid','yFgrid','yPdist','yTdist','yFdist','yPcumtrans',...
        'yPtrans','yPcumdist','yFcumdist','yTcumdist','ytrans',...
        'meannety','labtaxthresh','lumptransfer','ysortdist','ysort'};
    for i = 1:numel(newfields)
        income.(newfields{i}) = eval(newfields{i});
    end



end