function [avg_mpc1,avg_mpc4] = direct_mpcs(xgrid,p,income,model,prefs)
    % Computes MPCs by shifting the asset grid and comparing the new implied
    % consumption distribution to the stationary distribution of the new
    % asset grid

    ergodic_tol = 1e-11;
    ergodic_method = 2;
    
    if p.Display == 1
        disp('Using direct methods to find MPCs')
    end
    
    temp = cell(1,numel(p.mpcfrac));
    for im = 1:numel(p.mpcfrac)
        
        if p.mpcfrac{im} < 0 %|| p.mpcfrac{im} == 1e-5
            avg_mpc1(im) = NaN;
            avg_mpc2(im) = NaN;
            avg_mpc3(im) = NaN;
            avg_mpc4(im) = NaN;
            continue
        end
        
        % apply income transfer by adding small quantity to each point in 
        % asset space
        mpcamount       = p.mpcfrac{im} * income.meany * p.freq;
        mpcgrid_wide1   = xgrid.longgrid_wide + mpcamount;
        
        % add in new small xgrid values to replace some of those lost at
        % the bottom of grid
        newpts = 20;
        newvals = linspace(0,1,newpts+1)';
        newvals = newvals(1:end-1);
        
        % new grid
        mpc_grid_wide = zeros(p.nxlong+newpts,p.nyP,p.nyF);
        mpc_grid_wide(newpts+1:end,:,:) = mpcgrid_wide1;
        bottom = repmat(xgrid.longgrid_wide(1,:,:),[newpts 1 1]);
        mpc_grid_wide(1:newpts,:,:) = bottom + mpcamount * repmat(newvals,[1 p.nyP p.nyF]);
        
        % shift original stationary distribution by adding 0 probability for low x
        orig_SSdist                     = zeros(p.nxlong+newpts,p.nyP,p.nyF,p.nb);
        orig_SSdist(newpts+1:end,:,:,:) = model.SSdist_wide;
        
        % find stationary distribution of new asset space
        [mpcSSdist,mpctrans,SSsav,SScon] = find_stationary(p,model,income,...
                                        	prefs,mpc_grid_wide,ergodic_method,ergodic_tol);                        
                                            
        % compute average mpc's
        stationary_con = model.mean_c;
        avg_mpc1(im) = (orig_SSdist(:)'              * SScon(:) - stationary_con)/ mpcamount; 
        avg_mpc2(im) = (orig_SSdist(:)' * mpctrans   * SScon(:) - stationary_con)/ mpcamount + avg_mpc1(im);
        avg_mpc3(im) = (orig_SSdist(:)' * mpctrans^2 * SScon(:) - stationary_con)/ mpcamount + avg_mpc2(im);
        avg_mpc4(im) = (orig_SSdist(:)' * mpctrans^3 * SScon(:) - stationary_con)/ mpcamount + avg_mpc3(im);
    end

    

end