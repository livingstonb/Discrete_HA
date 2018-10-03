function [avg_mpc1,avg_mpc2] = direct_mpcs_new(xgrid,p,income,basemodel,prefs)

    NN = p.nxlong * p.nyP * p.nyF * p.nb;

    for im = 1:numel(p.mpcfrac)
        xx  = cell(1,4);
        sav = cell(1,4);
        con = cell(1,4);
        
        mpcamount = p.mpcfrac{1} * income.meany * p.freq;
        xx{1} = xgrid.longgrid_wide + mpcamount;
        
        %% One-period MPC
        
        con{1} = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            con{1}(:,iyP,iyF,ib) = basemodel.coninterp{iyP,iyF,ib}(xx{1}(:,iyP,iyF));
        end
        end
        end
        mpc{1} = (con{1}(:) - basemodel.con_longgrid) / mpcamount;
        avg_mpc{1} = basemodel.SSdist' * mpc{1};

        sav{1} = get_sav1(p,xx{1},basemodel);
        
        for t = 1:3
            % create matrix of xprime's, final dim NN by p.nyP*p.nyF*p.nyT
            % rows index points in today's asset space, columns index
            % next period's yP,yF,yT realizations
            xprime{t}   = zeros(NN,p.nyP*p.nyF*p.nb,p.nyT);
            for iyT = 1:p.nyT
                xprime{t}(:,:,iyT) = p.R * repmat(sav{t}(:),1,p.nyP*p.nyF)...
                                        + repmat(income.netymat(:,iyT)',NN,1);
            end
            xprime{t} = reshape(xprime{t},NN,p.nyP*p.nyF*p.nyT);

            expected_con = get_mpc(p,xprime{t},NN,basemodel,income,prefs,mpcamount);
            T{t} = transition_matrix(NN,xprime{t},p,xgrid,income,prefs);
            avg_mpc{2} = (basemodel.SSdist'* expected_con - basemodel.SSdist'*basemodel.con_longgrid) / mpcamount;
        end

        %% Multi-period MPCs
        for t = 2:4
            mpc{t} = get_mpc(p,x_mpc{t},NN,basemodel,income,prefs,mpcamount); 
            sav{t} = get_sav(p,x_mpc{t},basemodel,NN);
            x_mpc{t+1}  = find_xp(p,sav{t},income,NN);
        end
            
        
    end
    
   function sav = get_sav1(p,xx,basemodel)
        sav = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            sav(:,iyP,iyF,ib) = basemodel.savinterp{iyP,iyF,ib}(xx(:,iyP,iyF));
        end
        end
        end
   end

    function expected_con = get_mpc(p,x,NN,basemodel,income,prefs,mpcamount)
        con_mpc = zeros(NN,p.nyP*p.nyF*p.nb*p.nyT);
        x       = reshape(x,NN,p.nyP*p.nyF*p.nyT);
        rowblock = 1;
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            rowind1 = (rowblock-1)*p.nxlong + 1;
            rowind2 = rowblock*p.nxlong;
            x_rowblock = x(rowind1:rowind2,:);

            con_mpc_rowblock = basemodel.coninterp{iyP,iyF,ib}(x_rowblock(:));
            con_mpc(rowind1:rowind2,:) = reshape(con_mpc_rowblock,p.nxlong,p.nyP*p.nyF*p.nyT);

            rowblock = rowblock + 1;
        end
        end
        end

        % take expectation over yT
        ETcon_mpc = reshape(con_mpc,[],p.nyT) * income.yTdist;
        con_mpc  = reshape(ETcon_mpc,NN,p.nyP*p.nyF*p.nb);

        % construct transition matrix
        trans_exo = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
        trans_exo = kron(trans_exo,ones(p.nxlong,1));

        % expected consumption in each state (x,yP,yF,beta)
        expected_con = dot(trans_exo,con_mpc,2);

    end


    function T = transition_matrix(NN,xp,p,xgrid,income,prefs)
        T = zeros(NN,NN);
        col = 1;
        xp = reshape(xp,NN,p.nyP,p.nyF,p.nyT);
        for ib2 = 1:p.nb
        for iyF2 = 1:p.nyF
        for iyP2 = 1:p.nyP
            % create spline object
            fspace = fundef({'spli',xgrid.longgrid_wide(:,iyP2,iyF2),0,1});
            xp_T = reshape(squeeze(xp(:,iyP2,iyF2,:)),[],1);
            
            newcolumn = funbas(fspace,xp_T);
            newcolumn = reshape(newcolumn,[],p.nyT*p.nxlong) * kron(speye(p.nxlong),income.yTdist);
            T(:,p.nxlong*(col-1)+1:p.nxlong*col) = newcolumn;
            col = col + 1;
        end
        end
        end
        
        trans_exo = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
        trans_exo = kron(trans_exo,ones(p.nxlong));
        T = trans_exo .* T;
       
    end

end