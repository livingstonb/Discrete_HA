function [con,sav] = EGP_fun(muc,sgrid,xgrid,savtax,savtaxthreshold,...
                                u1inv,r)
% _wide variables have dimension nx by nyP*nyF*nb, or nx by N/nx

% consumption as a function of next period's cash (xprime)
con_s = u1inv(muc);
% cash-in-hand (x) as a function of xprime
x_s = sgrid + savtax * max(sgrid-savtaxthresh,0);

% interpolate from x(s) to get s(x)
x_s_wide = reshape(x_s,nx,N/nx);
sgrid_wide = reshape(sgrid,nx,N/nx);
V = interpOneD_vec(sgrid_wide,x_s_wide);

% saving as a function of x
sav_wide = V * sgrid_wide;

% deal with borrowing limit
constr_wide = xgrid < x_s_wide(1,:);
sav_wide(constr) = borrow_lim;
sav = reshape(sav_wide,N,1);

con = xgrid - sav;
end