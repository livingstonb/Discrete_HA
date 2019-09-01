classdef Income < handle
    properties (SetAccess = private)
        LoadIncome;
        Import;
        
        p;
        het;
        
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

        ytrans;
        meany1;

        ymatdist;
        labtaxthresh;

        lumptransfer;
        ymat;
        netymat;
        netymatEGP;
        netymatDST;
        meannety1;

        ysort;
        ysortdist;
        
        ytrans_live;
        ytrans_death;
    end
    
    methods
        function obj = Income(p,heterogeneity)
            obj.p = p;
            obj.het = heterogeneity;
            obj.LoadIncome = ~isempty(p.IncomeProcess);
            
            if obj.LoadIncome
                obj.Import = load(p.IncomeProcess);
            end
            
            obj.get_persistent_income();
            obj.get_transitory_income();
            obj.get_fixed_effect();
            obj.get_other_income_variables();
            
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
                    = aux.rouwenhorst(obj.p.nyP, -0.5*obj.p.sd_logyP^2, obj.p.sd_logyP, obj.p.rho_logyP);
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
                obj.logyTgrid = obj.Import.discmodel1.logyTgrid;
                obj.yTdist = obj.Import.discmodel1.yTdist;
                obj.p.nyT = length(obj.logyTgrid);
                obj.logyTgrid = reshape(obj.logyTgrid,[],1);
                obj.yTdist = reshape(obj.yTdist,[],1);
                obj.yTcumdist = cumsum(obj.yTdist,1);
            elseif obj.p.nyT>1
                %moments of mixture distribution
                lmu2 = obj.p.lambdaT.*obj.p.sd_logyT^2;
                lmu4 = 3.*obj.p.lambdaT.*(obj.p.sd_logyT^4);

                %fit those moments
                optionsNLLS = optimoptions('lsqnonlin','Display','Off');
                lpar = lsqnonlin(@(lp) aux.discretize_normal_var_kurt(...
                    lp,obj.p.nyT,-lmu2/2,lmu2,lmu4),[2 0.1],[],[],optionsNLLS);
                [lf,lx,lp] = aux.discretize_normal_var_kurt(lpar,obj.p.nyT,-lmu2/2,lmu2,lmu4);
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
                width = fzero(@(x) aux.discrete_normal(obj.p.nyF,-0.5*obj.p.sd_logyF^2 ,obj.p.sd_logyF ,x),2);
                [~,obj.logyFgrid,obj.yFdist] = aux.discrete_normal(obj.p.nyF,-0.5*obj.p.sd_logyF^2 ,obj.p.sd_logyF ,width);
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
            obj.ymat = repmat(obj.yPgrid,obj.p.nyF,1) ...
            	.* kron(obj.yFgrid,ones(obj.p.nyP,1)) * obj.yTgrid';

            % distribution of ymat
            obj.ymatdist = repmat(obj.yPdist,obj.p.nyF,1) ...
            	.* kron(obj.yFdist,ones(obj.p.nyP,1)) * obj.yTdist';

            % find mean y
            % isolate unique (yT,yF,yP) combinations
            temp = sortrows([obj.ymat(:) obj.ymatdist(:)],1);
            obj.ysort = temp(:,1);
            obj.ysortdist = temp(:,2);
            ycumdist_sort = cumsum(obj.ysortdist);

            % 1-period statistics
            obj.meany1 = obj.ymat(:)' * obj.ymatdist(:);
            totgrossy1 =obj. meany1;

            % find tax threshold on labor income
            if numel(obj.ysort)>1
                obj.labtaxthresh = lininterp1(ycumdist_sort,obj.ysort,obj.p.labtaxthreshpc);
            else
                obj.labtaxthresh = 0;
            end    

            % find net income
            totgrossyhigh = max(obj.ymat(:)-obj.labtaxthresh,0)' * obj.ymatdist(:);
            obj.lumptransfer = obj.p.labtaxlow * totgrossy1 ...
            	+ obj.p.labtaxhigh * totgrossyhigh;
            % netymat is N by nyT matrix
            obj.netymat = obj.lumptransfer + (1-obj.p.labtaxlow) * obj.ymat ...
            	- obj.p.labtaxhigh * max(obj.ymat-obj.labtaxthresh,0);
            obj.meannety1 = obj.netymat(:)' * obj.ymatdist(:);

            % net y values on HJB grid
            netymat_temp = reshape(obj.netymat,[1 obj.p.nyP obj.p.nyF obj.p.nyT]);
            obj.netymatEGP = repmat(netymat_temp,[obj.p.nx 1 1 1]);
            obj.netymatDST = repmat(netymat_temp,[obj.p.nx_DST 1 1 1]);

            % full transition matrix with beta and IES transitions, excluding and including death
            if obj.p.ResetIncomeUponDeath == 1
                yPtrans_death = repmat(obj.yPdist',obj.p.nyP,1);
            else
                yPtrans_death = obj.yPtrans;
            end

            if (numel(obj.p.risk_aver) == 1) && (numel(obj.p.invies) == 1) && (numel(obj.p.r)==1)
                obj.ytrans_live = kron(obj.het.betatrans,kron(eye(obj.p.nyF),obj.yPtrans));
                obj.ytrans_death = kron(obj.het.betatrans,kron(eye(obj.p.nyF),yPtrans_death));
            elseif numel(obj.p.r) > 1
                obj.ytrans_live = kron(obj.het.rtrans,kron(eye(obj.p.nyF),obj.yPtrans));
                obj.ytrans_death = kron(obj.het.rtrans,kron(eye(obj.p.nyF),yPtrans_death));
            else
                obj.ytrans_live = kron(obj.het.ztrans,kron(eye(obj.p.nyF),obj.yPtrans));
                obj.ytrans_death = kron(obj.het.ztrans,kron(eye(obj.p.nyF),yPtrans_death));
            end
        end
    end
end