function [con,sav] = find_policy(con,ns,nyP,nyF,nb,x_s,...
                            xgrid_wide,dieprob,betastacked,Emat,yTdist,...
                            sgrid_wide,savtax,savtaxthresh,...
                            u1inv,borrow_lim,N,max_iter,tol_iter)

iter = 1;
cdiff = 1;
% EGP ITERATION
while iter<max_iter && cdiff>tol_iter
    iter = iter + 1;
    conlast = con;

    % interpolate to get c(x') using c(x)
    conlast_wide = reshape(conlast,ns,nyP*nyF*nb);
    % initialize cons as function of x',yT
    c_xp = zeros(N,nyT);

    for iyT = 1:nyT
        x_s_wide = reshape(x_s(:,iyT),ns,nyP*nyF*nb); 
        c_xpT_wide = zeros(ns,nyP*nyF*nb);
        for col = 1:nyP*nyF*nb
            c_xpT_wide(:,col) = interp1(xgrid_wide(:,col),conlast_wide(:,col),x_s_wide(:,col),'linear','extrap');
        end
        c_xp(:,iyT)  = c_xpT_wide(:);
    end

    mucnext  = u1(c_xp);
    % muc this period as a function of s
    muc_s = (1-dieprob)*(1+r)*betastacked.*Emat*(mucnext*yTdist);
    % _wide variables have dimension nx by nyP*nyF*nb, or nx by N/nx

    % consumption as a function of next period's cash (xprime)
    con_s = u1inv(muc);
    % cash-in-hand (x) as a function of s
    x_s = sgrid_wide(:) + savtax * max(sgrid_wide(:)-savtaxthresh,0) + con_s;

    % interpolate from x(s) to get s(x), interpolate for each (beta,yP,yF)
    % separately
    x_s_wide = reshape(x_s,nx,N/nx);
    sav_wide = zeros(nx,N/nx);
    for col=1:N/nx
        sav_wide(:,col) = interp1(x_s_wide(:,col),sgrid_wide(:,col),xgrid_wide(:,col),'linear','extrap'); 
    end


    % deal with borrowing limit
    constr = xgrid_wide < x_s_wide(1,:);
    sav_wide(constr) = borrow_lim;
    sav = sav_wide(:);

    con = xgrid_wide(:) - sav - savtax * max(sav-savtaxthresh,0);

    cdiff = max(abs(con-conlast));
    if mod(iter,50) ==0
    	disp([' EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
    end
end

end