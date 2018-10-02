function [avg_mpc1,avg_mpc4] = direct_mpcs(xgrid,p,income,model,prefs)

    ergodic_tol = 1e-7;
    if p.Display == 1
        disp('Using direct methods to find MPCs')
    end
    
    model.a_longgrid
    xgrid.longgrid
    
    for im = 1:numel(p.mpcfrac)
        % apply income transfer by shifting grid
        mpcgrid = xgrid.longgrid + p.mpcfrac{im} * income.meany;
        
        % find stationary distribution associated with new grid
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            savlong(:,iyP,iyF,ib) = model.savinterp{iyP,iyF,ib}(xgrid.longgrid_wide(:,iyP,iyF));
        end
        end
        end
        
        % use stationary function
    end

    

end