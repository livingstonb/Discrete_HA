classdef Income < handle
    properties (SetAccess = private)
        LoadIncome;
        Import;
        
        p;
        
        logyPgrid;
        yPdist;
        yPtrans;
        yPgrid
        yPcumdist;
        yPcumtrans;
        
        logyTgrid;
        yTdist;
        yTgrid;
        yTcumdist;
        
        logyFgrid;
        yFgrid;
        yFdist;
        yFcumdist;
    end
    
    methods
        function obj = Income(p,heterogeneity)
            obj.p = p;
            obj.LoadIncome = ~isempty(p.IncomeProcess);
            
            if obj.LoadIncome
                obj.Import = load(p.IncomeProcess);
            end
            
            obj.get_persistent_income();
            obj.get_transitory_income();
            
            if size(obj.yTgrid,2)>1 || size(obj.yFgrid,2)>1 || size(obj.yPgrid,2)>1
                error('All income grids must be column vectors')
            end
            if size(obj.yTdist,2)>1 || size(obj.yFdist,2)>1 || size(obj.yPdist,2)>1
                error('All income distributions must be column vectors')
            end
        end
        
        function get_persistent_income(obj)
            if obj.LoadIncome
                obj.logyPgrid = obj.Import.discmodel1.logyPgrid;
                obj.yPdist = obj.Import.discmodel1.yPdist;
                obj.yPtrans = obj.Import.discmodel1.yPtrans;
                obj.p.nyP = length(obj.logyPgrid);
                obj.logyPgrid = reshape(obj.logyPgrid,[],1);
                obj.yPdist = reshape(obj.yPdist,[],1);
            elseif obj.p.nyP > 1
                [obj.logyPgrid, obj.yPtrans, obj.yPdist] ...
                    = rouwenhorst(obj.p.nyP, -0.5*obj.p.sd_logyP^2, obj.p.sd_logyP, obj.p.rho_logyP);
            else
                obj.logyPgrid = 0;
                obj.yPdist = 1;
                obj.yPtrans = 1;
            end
            
            obj.yPgrid = exp(obj.logyPgrid);
            obj.yPgrid = obj.yPgrid/(obj.yPdist'*obj.yPgrid);
            obj.logyPgrid = log(obj.yPgrid);
            obj.yPcumdist = cumsum(obj.yPdist,1);
            obj.yPcumtrans = cumsum(obj.yPtrans,2);
        end
        
        function get_transitory_income(obj)
            if obj.LoadIncome
                obj.logyTgrid = Import.discmodel1.logyTgrid;
                obj.yTdist = Import.discmodel1.yTdist;
                obj.p.nyT = length(obj.logyTgrid);
                obj.logyTgrid = reshape(obj.logyTgrid,[],1);
                obj.yTdist = reshape(obj.yTdist,[],1);
                obj.yTcumdist = cumsum(obj.yTdist,1);
            elseif obj.p.nyT>1
                %moments of mixture distribution
                lmu2 = obj.p.lambdaT.*obj.p.sd_logyT^2;
                lmu4 = 3.*obj.p.lambdaT.*(obj.p.sd_logyT^4);

                %fit thjose moments
                optionsNLLS = optimoptions(@lsqnonlin,'Display','Off');
                lpar = lsqnonlin(@(lp)discretize_normal_var_kurt(...
                    lp,obj.p.nyT,-lmu2/2,lmu2,lmu4),[2 0.1],[],[],optionsNLLS);
                [lf,lx,lp] = discretize_normal_var_kurt(lpar,obj.p.nyT,-lmu2/2,lmu2,lmu4);
                obj.logyTgrid = lx;
                obj.yTdist = lp;
                obj.yTcumdist = cumsum(obj.yTdist,1);

            elseif obj.p.nyT==1
                obj.logyTgrid = 0;
                obj.yTdist = 1;
                obj.yTcumdist = 1;
            end
            
            obj.yTgrid = exp(obj.logyTgrid);
            obj.yTgrid = obj.yTgrid / (obj.yTdist'*obj.yTgrid);
            obj.logyTgrid = log(obj.yTgrid);
        end
        
        function get_fixed_effect(obj)
            if obj.p.nyF>1
                width = fzero(@(x)discrete_normal(obj.p.nyF,-0.5*obj.p.sd_logyF^2 ,obj.p.sd_logyF ,x),2);
                [~,obj.logyFgrid,obj.yFdist] = discrete_normal(obj.p.nyF,-0.5*obj.p.sd_logyF^2 ,obj.p.sd_logyF ,width);
                obj.logyFgrid = reshape(obj.logyFgrid,[],1);
                obj.yFdist = reshape(obj.yFdist,[],1);
            elseif obj.p.nyF==1
                obj.logyFgrid = 0;
                obj.yFdist = 1;
            end
            obj.yFgrid = exp(obj.logyFgrid);
            % normalize fixed effect such that mean = 1 if annual, 1/4 if quarterly
            obj.yFgrid = obj.yFgrid/(obj.yFdist'*obj.yFgrid*obj.p.freq);
            obj.logyFgrid = log(obj.yFgrid);
            obj.yFcumdist = cumsum(obj.yFdist,1);
        end
        
        function get_other_income_variables(obj)
            % transition probabilities for yP-yF combined grid
            obj.ytrans = kron(eye(obj.p.nyF),obj.yPtrans);

            % construct matrix of y combinations
            obj.ymat = repmat(yPgrid,p.nyF,1) .* kron(yFgrid,ones(p.nyP,1)) * yTgrid';

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

            if (numel(p.risk_aver) == 1) && (numel(p.invies) == 1) && (numel(p.r)==1)
                ytrans_live = kron(heterogeneity.betatrans,kron(eye(p.nyF),yPtrans));
                ytrans_death = kron(heterogeneity.betatrans,kron(eye(p.nyF),yPtrans_death));
            elseif numel(p.r) > 1
                ytrans_live = kron(heterogeneity.rtrans,kron(eye(p.nyF),yPtrans));
                ytrans_death = kron(heterogeneity.rtrans,kron(eye(p.nyF),yPtrans_death));
            else
                ytrans_live = kron(heterogeneity.ztrans,kron(eye(p.nyF),yPtrans));
                ytrans_death = kron(heterogeneity.ztrans,kron(eye(p.nyF),yPtrans_death));
            end
        end
    end
end