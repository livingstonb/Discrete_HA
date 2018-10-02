function [avg_mpc1,avg_mpc4] = direct_mpcs(xgrid,p,income,model,prefs)

    ergodic_tol = 1e-7;
    if p.Display == 1
        disp('Using direct methods to find MPCs')
    end
    
    for im = 1:numel(p.mpcfrac)
        
        if p.mpcfrac{im} < 0
            avg_mpc1(im) = NaN;
            avg_mpc2(im) = NaN;
            avg_mpc3(im) = NaN;
            avg_mpc4(im) = NaN;
            continue
        end
        
        % apply income transfer by shifting grid
        mpcamount       = p.mpcfrac{im} * income.meany;
        mpcgrid_wide1   = xgrid.longgrid_wide + mpcamount;
        % add in new small xgrid values to replace some of those missing
        newvals = linspace(0,1,21)';
        newvals = newvals .^ (1/0.5);
        newvals = newvals(1:end-1);
        
        mpc_grid_wide = zeros(p.nxlong+20,p.nyP,p.nyF);
        mpc_grid_wide(21:end,:,:) = mpcgrid_wide1;
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            bottom  = xgrid.longgrid_wide(1,iyP,iyF);
            mpc_grid_wide(1:20,iyP,iyF) = bottom + mpcamount * newvals;
        end
        end
        
        % shift original stationary distribution by adding 0 probability for low x
        orig_SSdist                 = zeros(p.nxlong+20,p.nyP,p.nyF,p.nb);
        orig_SSdist(21:end,:,:,:)   = model.SSdist_wide;
        
        % find stationary distribution of new asset space
        [mpcSSdist,mpctrans,SSsav,SScon] = find_stationary(p,model,income,...
                                                prefs,mpc_grid_wide,ergodic_tol);
                                            
        
        avg_mpc1(im) = (orig_SSdist(:)' - mpcSSdist') * SScon(:) / mpcamount; 
        avg_mpc2(im) = (orig_SSdist(:)'*mpctrans   - mpcSSdist') * SScon(:) / mpcamount + avg_mpc1(im);
        avg_mpc3(im) = (orig_SSdist(:)'*mpctrans^2 - mpcSSdist') * SScon(:) / mpcamount + avg_mpc2(im);
        avg_mpc4(im) = (orig_SSdist(:)'*mpctrans^3 - mpcSSdist') * SScon(:) / mpcamount + avg_mpc3(im);
    end

    

end