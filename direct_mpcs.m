function [avg_mpc1,avg_mpc4] = direct_mpcs(xgrid,p,income,model,prefs)

    ergodic_tol = 1e-7;
    if p.Display == 1
        disp('Using direct methods to find MPCs')
    end
    
    for im = 1:numel(p.mpcfrac)
        %% GET STATIONARY DISTRIBUTION OF SHIFTED GRID
        
        % apply income transfer by shifting grid
        mpcgrid = xgrid.longgrid_wide + p.mpcfrac{im} * income.meany;
        
        % add xgrid minimums back to new grid
        % and shift original stationary distribution by adding 0 prob to
        % minimums of xgrid
        xgrid_mins              = min(xgrid.longgrid_wide,[],1);
        mpcgridnew              = zeros(p.nxlong+1,p.nyP,p.nyF,p.nb);
        mpcgridnew(2:end,:,:) = mpcgrid;
        mpcgrid = mpcgridnew;
        mpcgrid(1,:,:) = xgrid_mins;
        % shift original stationary distribution by adding 0 for xgrid mins
        orig_SS_dist            = zeros(p.nxlong+1,p.nyP,p.nyF,p.nb);
        orig_SS_dist(2:end,:,:,:) = model.SS_dist_wide;
%         for iyF = 1:p.nyF
%         for iyP = 1:p.nyP
%             mpcgrid(1,iyF,iyP) = xgrid_mins(iyF,iyP);
%             
%         end
%         end

        % find policy function interpolated onto new grid
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            sav(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(mpcgrid(:,iyP,iyF));
            con(:,iyP,iyF,ib) = model.coninterp{iyP,iyF,ib}(mpcgrid(:,iyP,iyF));
        end
        end
        end
        
        % find stationary distribution of new asset space
        [mpcSSdist,mpctrans,SSsav,SScon] = find_stationary(p,model,income,...
                                    prefs,xgrid.longgrid_wide,ergodic_tol);
                                
        mpc1{im} = (orig_SS_dist - mpcSSdist)' * 
    end

    

end