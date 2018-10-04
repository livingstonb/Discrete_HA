function [beta,exitflag] = iterate_beta(p,xgrid,sgrid,prefs,income)
    % This function is used to speed up convergence to the correct beta.
    % fzero first optimizes using solve_EGP with a small grid size and high
    % value for ergodic_tol. After convergence, this function shrinks the search 
    % interval for beta and calls fzero while passing solve_EGP a larger 
    % grid size and a smaller value for ergodic_tol.
    
    
    %% LOW TOLERANCE, SMALL GRID
    
    % Set interval for iteration
    if p.nb == 1
        beta_ub = p.betaH - 1e-5;
    else
        beta_ub = p.betaH - 1e-5 - betawidth;
    end
    beta_lb = p.betaL;

    if p.FastIter == 1
        % Narrow the interval by optimizing over a small grid with low
        % tolerance
        ergodic_tol = 1e-5;
        ergodic_method = 1;
        Intermediate = 1;
        iterate_EGP = @(x) solve_EGP(x,p,...
                xgrid,sgrid,prefs,ergodic_method,ergodic_tol,income,Intermediate);
    
        low_tol = 1e-4;
        check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
        options = optimset('TolX',low_tol,'OutputFcn',check_evals);
        [beta,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);

        if exitflag ~= 1
            return
        end
        
        beta_lb = beta - 0.01;
        beta_ub = beta + 0.01;
    end
        
    %% HIGH TOLERANCE, LARGE GRID
    fprintf('Switching to higher tolerance\n')
    ergodic_tol = 1e-7;
    ergodic_method = 1;
    Intermediate = 0;
    iterate_EGP = @(x) solve_EGP(x,p,...
        xgrid,sgrid,prefs,ergodic_method,ergodic_tol,income,Intermediate);

    high_tol = p.tolAY;
    check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
    options = optimset('TolX',high_tol,'OutputFcn',check_evals);
    [beta,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);
    end
