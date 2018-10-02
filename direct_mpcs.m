function [avg_mpc1,avg_mpc4] = direct_mpcs(xgrid,p,income,model,prefs)

    ergodic_tol = 1e-7;
    if p.Display == 1
        disp('Using direct methods to find MPCs')
    end
    
    for im = 1:numel(p.mpcfrac)
        
        if p.mpcfrac{im} < 0
            avg_mpc1(im) = 0;
            avg_mpc2(im) = 0;
            continue
        end
        
        % apply income transfer by shifting grid
        xgrid_mins      = min(xgrid.longgrid_wide,[],1);
        mpcamount       = p.mpcfrac{im} * income.meany;
        mpcgrid_wide    = xgrid.longgrid_wide + mpcamount;
        
        % add xgrid minimums back to new grid
        % and shift original stationary distribution by adding 0 prob to
        % minimums of xgrid
        mpcgridnew              = zeros(p.nxlong+1,p.nyP,p.nyF,p.nb);
        mpcgridnew(2:end,:,:)   = mpcgrid_wide;
        mpcgrid_wide            = mpcgridnew;
        mpcgrid_wide(1,:,:)     = xgrid_mins;
        
        % shift original stationary distribution by adding 0 for xgrid mins
        orig_SSdist                = zeros(p.nxlong+1,p.nyP,p.nyF,p.nb);
        orig_SSdist(2:end,:,:,:)   = model.SSdist_wide;
        
        % find stationary distribution of new asset space
        [mpcSSdist,mpctrans,SSsav,SScon] = find_stationary(p,model,income,...
                                                prefs,mpcgrid_wide,ergodic_tol);
                        
        avg_mpc1(im) = (orig_SSdist(:)' - mpcSSdist') * SScon(:) / mpcamount; 
        avg_mpc2(im) = (orig_SSdist(:)'*mpctrans   - mpcSSdist') * SScon(:) / mpcamount + avg_mpc1(im);
        avg_mpc3(im) = (orig_SSdist(:)'*mpctrans^2 - mpcSSdist') * SScon(:) / mpcamount + avg_mpc2(im);
        avg_mpc4(im) = (orig_SSdist(:)'*mpctrans^3 - mpcSSdist') * SScon(:) / mpcamount + avg_mpc3(im);
    end

    

end