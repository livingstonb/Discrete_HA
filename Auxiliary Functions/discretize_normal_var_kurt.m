function [f,lx,lp] = discretize_normal_var_kurt(y,n,mu2,mu4)
    lwidth = y(1);
    llambda = y(2);
    lx = linspace(-lwidth*sqrt(mu2), lwidth*sqrt(mu2),n)';
    lp = discrete_normal_alt(lx,0,sqrt(mu2));
    lmass0 = zeros(n,1);
    lmass0((n+1)/2) = 1; 
    
    lp = (1-llambda).*lp + lmass0.*llambda;
    
    Ex2 = sum(lp.*(lx.^2));
    Ex4 = sum(lp.*(lx.^4));
    
    f = [sqrt(Ex2)-sqrt(mu2); Ex4.^0.25 - mu4.^0.25];
    
end 

