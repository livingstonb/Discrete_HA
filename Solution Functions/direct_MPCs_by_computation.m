function [MPCs,agrid_dist,norisk_mpcs1_a_direct]...
                = direct_MPCs_by_computation(p,basemodel,income,prefs,agrid_short,norisk)

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
    
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im)*income.meany1*p.freq;
        end

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

        if im == 0
            con_baseline = Econ;
        else
            if p.Display == 1
                disp('Computing IMPC(1,1)')
            end
            % Compute m(a,yP,yF,beta) = E[m(x,yP,yF,beta)|a,yP,yF,beta]
            % MPC in period 1 out of period 1 shock
            mpcs_1_1_yP_yF_beta = (Econ - con_baseline) / mpcamount;
            MPCs.avg_1_1(im) = basemodel.adist(:)' * mpcs_1_1_yP_yF_beta(:);

            % Compute m(a) = E(m(a,yP,yF,beta)|a)
            %       = sum of P(yP,yF,beta|a) * m(a,yP,yF,beta) over all
            %         (yP,yF,beta)
            MPCs.mpcs_1_1{im} = sum(sum(sum(Pcondl .* mpcs_1_1_yP_yF_beta,4),3),2);
        end
        
        % 4-PERIOD MPCS
        % Construct transition matrix from period 1 to period 2
        NN = p.nxlong*p.nyP*p.nyF*p.nb;
        if im>0

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

            % MPC in period 2 out of period 1 shock
            if p.Display == 1
                disp('Computing IMPC(1,2)')
            end
            mpcs_1_2_yP_yF_beta = (T12-speye(NN))*con_baseline(:)/mpcamount;
            mpcs_1_2_yP_yF_beta = reshape(mpcs_1_2_yP_yF_beta,[p.nxlong p.nyP p.nyF p.nb]);
            MPCs.mpcs_1_2{im} = sum(sum(sum(Pcondl .* mpcs_1_2_yP_yF_beta,4),3),2);
            MPCs.avg_1_2(im) = basemodel.adist(:)' * mpcs_1_2_yP_yF_beta(:);

            % MPC in period 3 out of period 1 shock
            if p.Display == 1
                disp('Computing IMPC(1,3)')
            end
            mpcs_1_3_yP_yF_beta = (T12*basemodel.statetrans-speye(NN))*con_baseline(:)/mpcamount;
            mpcs_1_3_yP_yF_beta = reshape(mpcs_1_3_yP_yF_beta,[p.nxlong p.nyP p.nyF p.nb]);
            MPCs.mpcs_1_3{im} = sum(sum(sum(Pcondl .* mpcs_1_3_yP_yF_beta,4),3),2);
            MPCs.avg_1_3(im) = basemodel.adist(:)' * mpcs_1_3_yP_yF_beta(:);

            % MPC in period 4 out of period 1 shock
            if p.Display == 1
                disp('Computing IMPC(1,4)')
            end
            mpcs_1_4_yP_yF_beta = (T12*basemodel.statetrans^2-speye(NN))*con_baseline(:)/mpcamount;
            mpcs_1_4_yP_yF_beta = reshape(mpcs_1_4_yP_yF_beta,[p.nxlong p.nyP p.nyF p.nb]);
            MPCs.mpcs_1_4{im} = sum(sum(sum(Pcondl .* mpcs_1_4_yP_yF_beta,4),3),2);
            MPCs.avg_1_4(im) = basemodel.adist(:)' * mpcs_1_4_yP_yF_beta(:);
            
            if p.freq==1
                MPCs.avg_1_1to4(im) = NaN;
                MPCs.avg_1_5to8(im) = NaN;
                MPCs.avg_1_9to12(im) = NaN;
                MPCs.avg_1_13to16(im) = NaN;
            else
                % MPC in year 1 out of quarter 1 shock
                if p.Display == 1
                    disp('Computing IMPC(1,1-4)')
                end
                MPCs.avg_1_1to4(im) = MPCs.avg_1_1(im) + MPCs.avg_1_2(im) + MPCs.avg_1_3(im) + MPCs.avg_1_4(im);
                
                % MPC in periods 5-8 out of period 1 shock
                if p.Display == 1
                    disp('Computing IMPC(1,5-8)')
                end
                mpcs_1_5_yP_yF_beta = (T12*basemodel.statetrans^3-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_6_yP_yF_beta = (T12*basemodel.statetrans^4-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_7_yP_yF_beta = (T12*basemodel.statetrans^5-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_8_yP_yF_beta = (T12*basemodel.statetrans^5-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_5to8_yP_yF_beta = mpcs_1_5_yP_yF_beta+mpcs_1_6_yP_yF_beta+...
                                        mpcs_1_7_yP_yF_beta+mpcs_1_8_yP_yF_beta;
                MPCs.avg_1_5to8(im) = basemodel.adist(:)' * mpcs_1_5to8_yP_yF_beta(:);
                clear mpcs_1_5_yP_yF_beta mpcs_1_6_yP_yF_beta mpcs_1_7_yP_yF_beta
                clear mpcs_1_8_yP_yF_beta mpcs_1_5to8_yP_yF_beta

                % MPC in periods 9-12 out of period 1 shock
                if p.Display == 1
                    disp('Computing IMPC(1,9-12)')
                end
                mpcs_1_9_yP_yF_beta = (T12*basemodel.statetrans^6-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_10_yP_yF_beta = (T12*basemodel.statetrans^7-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_11_yP_yF_beta = (T12*basemodel.statetrans^8-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_12_yP_yF_beta = (T12*basemodel.statetrans^9-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_9to12_yP_yF_beta = mpcs_1_9_yP_yF_beta+mpcs_1_10_yP_yF_beta+...
                                        mpcs_1_11_yP_yF_beta+mpcs_1_12_yP_yF_beta;
                MPCs.avg_1_9to12(im) = basemodel.adist(:)' * mpcs_1_9to12_yP_yF_beta(:);
                clear mpcs_1_9_yP_yF_beta mpcs_1_10_yP_yF_beta mpcs_1_11_yP_yF_beta...
                        mpcs_1_12_yP_yF_beta mpcs_1_9to12_yP_yF_beta

                % MPC in periods 13-16 out of period 1 shock
                if p.Display == 1
                    disp('Computing IMPC(1,13-16)')
                end
                mpcs_1_13_yP_yF_beta = (T12*basemodel.statetrans^10-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_14_yP_yF_beta = (T12*basemodel.statetrans^11-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_15_yP_yF_beta = (T12*basemodel.statetrans^12-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_16_yP_yF_beta = (T12*basemodel.statetrans^13-speye(NN))*con_baseline(:)/mpcamount;
                mpcs_1_13to16_yP_yF_beta = mpcs_1_13_yP_yF_beta+mpcs_1_14_yP_yF_beta+...
                                        mpcs_1_15_yP_yF_beta+mpcs_1_16_yP_yF_beta;
                MPCs.avg_1_13to16(im) = basemodel.adist(:)' * mpcs_1_13to16_yP_yF_beta(:);
                clear mpcs_1_13_yP_yF_beta mpcs_1_14_yP_yF_beta mpcs_1_15_yP_yF_beta...
                        mpcs_1_16_yP_yF_beta mpcs_1_13to16_yP_yF_beta
            end
        end
    end

    % Distribution over agrid, P(a)
    agrid_dist = sum(sum(sum(basemodel.adist,4),3),2);
    
    %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITHOUT INCOME RISK)
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im)*income.meany1*p.freq;
        end
        
        x_mpc = agrid_short + income.meany1 + mpcamount;
        con = zeros(p.nxlong,p.nb);
        for ib = 1:p.nb
            con(:,ib) = norisk.coninterp{ib}(x_mpc);
        end
        
        if im == 0
            con_baseline = con;
        else
            % Compute m(a,beta)
            mpcs1_a_beta = (con - con_baseline) / mpcamount;

            % Compute m(x) = E(m(x,beta)|x)
            %       = sum of P(beta|x) * m(x,beta) over all beta
            % beta is exogenous so P(beta|x) = P(beta)
            norisk_mpcs1_a_direct{im} = mpcs1_a_beta * prefs.betadist;
        end
    end
end