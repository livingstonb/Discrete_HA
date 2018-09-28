function [beta,exitflag] = iterate_beta(p,xgrid,sgrid,prefs,income)
    %% LOW TOLERANCE, SMALL GRID
    ergodic_tol = 1e-4;
    iterate_EGP = @(x) solve_EGP(x,p,...
    xgrid,sgrid,prefs,ergodic_tol,income,50);

    beta_lb = p.betaL;
    if p.nb == 1
        beta_ub = p.betaH - 1e-5;
    else
        beta_ub = p.betaH - 1e-5 - betawidth;
    end

    low_tol = 1e-3;
    check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
    options = optimset('TolX',low_tol,'OutputFcn',check_evals);
    [beta,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);

    if exitflag ~= 1
        return
    end
        
    %% HIGH TOLERANCE, LARGER GRID
    fprintf('Switching to higher tolerance\n')
    ergodic_tol = 1e-6;
    ExpandGrid = 1;
    iterate_EGP = @(x) solve_EGP(x,p,...
    xgrid,sgrid,prefs,ergodic_tol,income,p.nxlong_int);

    beta_lb = beta - 0.01;
    beta_ub = beta + 0.01;

    high_tol = 1e-5;
    check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
    options = optimset('TolX',high_tol,'OutputFcn',check_evals);
    [beta,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);
    end
