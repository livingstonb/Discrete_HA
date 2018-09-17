function [con,sav] = EGP_fun(muc,sgrid,sgrid_long,xgrid,savtax,savtaxthresh,...
                                u1inv,r,borrow_lim,N)
% _wide variables have dimension nx by nyP*nyF*nb, or nx by N/nx
nx = size(xgrid,1);

% consumption as a function of next period's cash (xprime)
con_s = u1inv(muc);
% cash-in-hand (x) as a function of xprime
x_s = sgrid_long + savtax * max(sgrid_long-savtaxthresh,0) + con_s;

% interpolate from x(s) to get s(x)
x_s_wide = reshape(x_s,nx,N/nx);
sgrid_wide = repmat(sgrid,1,N/nx);
save_wide = zeros(nx,N/nx);
for col=1:N/nx
    V = interpOneD_vec(sgrid_wide(:,col),x_s_wide(:,col));
    save_wide(:,col) = V * sgrid_wide(:,col);
end

% saving as a function of x
sav_wide = V * sgrid_wide;

% deal with borrowing limit
constr = xgrid < x_s_wide(1,:);
sav_wide(constr) = borrow_lim;
sav = reshape(sav_wide,N,1);

con = repmat(xgrid,N/nx,1) - sav;
end