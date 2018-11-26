function [MPCs,agrid_dist] = direct_MPCs_s_t(p,basemodel,models,income,prefs,agrid_short)

    if p.Display == 1
        disp('Computing MPCs')
    end
    
    %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITH INCOME RISK)

    % Find P(yP,yF,beta|a) = P(a,yP,yF,beta)/P(a)
    Pa = sum(sum(sum(basemodel.adist,4),3),2);
    Pa = repmat(Pa,[1 p.nyP p.nyF p.nb]);
    Pcondl = basemodel.adist ./ Pa;
    Pcondl(Pa == 0) = 0;

    % Each (a,yP,yF) is associated with nyT possible x values, create this
    % grid here
    netymat_reshape = reshape(income.netymat,[1 p.nyP p.nyF p.nyT]);
    netymat_reshape = repmat(netymat_reshape,[p.nxlong 1 1 1]);
    xgrid_yT = repmat(agrid_short,[1 p.nyP p.nyF p.nyT]) + netymat_reshape;

    % For 4-period MPCs
    fspace = fundef({'spli',agrid_short,0,1});
    % transition matrix between (yP,yF,beta) states, cond'l on living
    trans_live = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    % transition matrix between (yP,yF,beta) states cond'l on dying
    yPtrans_stationary = repmat(income.yPdist',p.nyP,1);
    trans_death = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans_stationary));
    
    for mpcamount = [0, 0.01*income.meany1*p.freq]

        x_mpc = xgrid_yT + mpcamount;
        sav = zeros(p.nxlong,p.nyP,p.nyF,p.nb,p.nyT);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,:);
            sav_iyP_iyF_iyT = basemodel.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
            sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,[p.nxlong 1 1 1 p.nyT]);
        end
        end
        end
        sav = max(sav,p.borrow_lim);

        % Take expectation over yT
        % con becomes E[con(x,yP,yF,beta)|a,yP,yF,beta]
        Esav = reshape(sav,[],p.nyT) * income.yTdist;
        Esav = reshape(Esav,[p.nxlong p.nyP p.nyF p.nb]);
        Exx = reshape(reshape(x_mpc,[],p.nyT) * income.yTdist,[p.nxlong p.nyP p.nyF]);
        Econ = repmat(Exx,[1 1 1 p.nb]) - Esav - p.savtax * max(Esav - p.savtaxthresh,0);  
        clear Esav Exx 

        if mpcamount == 0
            con_baseline = Econ;
        else
            if p.Display == 1
                disp('Computing IMPC(1,1)')
            end
            % Compute m(a,yP,yF,beta) = E[m(x,yP,yF,beta)|a,yP,yF,beta]
            % MPC in period 1 out of period 1 shock
            mpcs_1_1_yP_yF_beta = (Econ - con_baseline) / mpcamount;
            MPCs.avg_s_t{1,1} = basemodel.adist(:)' * mpcs_1_1_yP_yF_beta(:);

            % Compute m(a) = E(m(a,yP,yF,beta)|a)
            %       = sum of P(yP,yF,beta|a) * m(a,yP,yF,beta) over all
            %         (yP,yF,beta)
%             mpcs_s_t{1,1} = sum(sum(sum(Pcondl .* mpcs_1_1_yP_yF_beta,4),3),2);
        end
        
        % 4-PERIOD AND ABOVE MPCS
        % Construct transition matrix from period 1 to period 2
        
        if mpcamount>0
            NN = p.nxlong*p.nyP*p.nyF*p.nb;
            % transition from receiving MPC back to regime of usual policy function (T12)
            aprime_live = p.R * sav;
            interp = funbas(fspace,aprime_live(:));
            interp = reshape(interp,NN,p.nxlong*p.nyT) * kron(speye(p.nxlong),income.yTdist);
            if p.Bequests == 1
                interp_death = interp;
            else
                interp_death = sparse(NN,p.nxlong);
                interp_death(:,1) = 1;
            end
            
            T12 = sparse(NN,NN);
            for col = 1:p.nyP*p.nyF*p.nb
                newblock_live = bsxfun(@times,kron(trans_live(:,col),ones(p.nxlong,1)),interp);
                newblock_death = bsxfun(@times,kron(trans_death(:,col),ones(p.nxlong,1)),interp_death);
                T12(:,p.nxlong*(col-1)+1:p.nxlong*col) = (1-p.dieprob)*newblock_live + p.dieprob*newblock_death;
            end
            clear interp_death interp_live newblock_death newblock_death

            % FIND CONSUMPTION 
            for is = 2:4
                for it = 1:is-1 % for t < s for now                
                    for ii = 1:it
                        % construct transition matrix from 'ii' to 'ii' + 1
                        x_mpc = xgrid_yT;
                        sav = zeros(p.nxlong,p.nyP,p.nyF,p.nb,p.nyT);
                        for ib = 1:p.nb
                        for iyF = 1:p.nyF
                        for iyP = 1:p.nyP
                            x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,:);
                            sav_iyP_iyF_iyT = models{is,ii}.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
                            sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,[p.nxlong 1 1 1 p.nyT]);
                        end
                        end
                        end
                        sav = max(sav,p.borrow_lim);

                        aprime_live = p.R * sav;
                        interp = funbas(fspace,aprime_live(:));
                        interp = reshape(interp,NN,p.nxlong*p.nyT) * kron(speye(p.nxlong),income.yTdist);
                        if p.Bequests == 1
                            interp_death = interp;
                        else
                            interp_death = sparse(NN,p.nxlong);
                            interp_death(:,1) = 1;
                        end

                        % Transition matrix from 'ii' to 'ii'+1
                        Ti = sparse(NN,NN);
                        for col = 1:p.nyP*p.nyF*p.nb
                            newblock_live = bsxfun(@times,kron(trans_live(:,col),ones(p.nxlong,1)),interp);
                            newblock_death = bsxfun(@times,kron(trans_death(:,col),ones(p.nxlong,1)),interp_death);
                            Ti(:,p.nxlong*(col-1)+1:p.nxlong*col) = (1-p.dieprob)*newblock_live + p.dieprob*newblock_death;
                        end
                        clear interp_death interp_live newblock_death newblock_live

                        % Create transition matrix from period 1 to period t
                        if ii == 1
                            T1t = speye(NN);
                        else
                            T1t = T1t * Ti;
                        end
                    end
                    
                    % Take expectation over yT
                    Esav = reshape(sav,[],p.nyT) * income.yTdist;
                    Esav = reshape(Esav,[p.nxlong p.nyP p.nyF p.nb]);
                    Exx = reshape(reshape(x_mpc,[],p.nyT) * income.yTdist,[p.nxlong p.nyP p.nyF]);
                    Econ = repmat(Exx,[1 1 1 p.nb]) - Esav - p.savtax * max(Esav - p.savtaxthresh,0);  
                    clear Esav Exx 
                    con = Econ;

                    % now compute IMPC(s,t)
                    mpcs = ( T1t * con(:) - con_baseline(:) ) / mpcamount;
                    MPCs.avg_s_t{is,it} = basemodel.adist(:)' * mpcs(:);
                    
                end
                
                it = is;
            end
        end
    end



    % Distribution over agrid, P(a)
    agrid_dist = sum(sum(sum(basemodel.adist,4),3),2);
    
   
end