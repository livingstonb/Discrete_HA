function [con,sav] = EGP_fun(muc,sgrid_wide,xgrid_wide,savtax,savtaxthresh,...
                                u1inv,borrow_lim,N)
% _wide variables have dimension nx by nyP*nyF*nb, or nx by N/nx
nx = size(xgrid_wide,1);

% consumption as a function of next period's cash (xprime)
con_s = u1inv(muc);
% cash-in-hand (x) as a function of s
x_s = sgrid_wide(:) + savtax * max(sgrid_wide(:)-savtaxthresh,0) + con_s;

% interpolate from x(s) to get s(x), interpolate for each (beta,yP,yF)
% separately
x_s_wide = reshape(x_s,nx,N/nx);
sav_wide = zeros(nx,N/nx);
for col=1:N/nx
    %V = interpOneD_vec(sgrid_wide(:,col),x_s_wide(:,col));
    %sav_wide(:,col) = V * sgrid_wide(:,col);
    sav_wide(:,col) = lininterp1(x_s_wide(:,col),sgrid_wide(:,col),xgrid_wide(:,col)); 
end


% deal with borrowing limit
constr = xgrid_wide < x_s_wide(1,:);
sav_wide(constr) = borrow_lim;
sav = sav_wide(:);

con = xgrid_wide(:) - sav - savtax * max(sav-savtaxthresh,0);
end