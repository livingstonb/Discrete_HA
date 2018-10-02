function makeplots(p,xgrid,sgrid,basemodel,income,sim_results,assetmeans)
    % This function makes plots based on fields from 'basemodel' 'income'
    % and 'simulations'

    figure(1);

    %plot for median fixed effect
    if mod(p.nyF,2)==1
        iyF = (p.nyF+1)/2;
    else
        iyF = p.nyF/2;
    end

    % plot for first beta
    iyb = 1;
    % if nb = 1, force plot of first beta
    if p.nb == 1
        iyb = 1;
    end

    % consumption policy function
    subplot(2,4,1);
    plot(xgrid.orig_wide(:,1,iyF,iyb),basemodel.con_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.con_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',1);
    grid;
    xlim([p.borrow_lim p.xmax]);
    title('Consumption Policy Function');
    legend('Lowest income state','Highest income state');

    % savings policy function
    subplot(2,4,2);
    plot(xgrid.orig_wide(:,1,iyF,iyb),basemodel.sav_wide(:,1,iyF,iyb)./xgrid.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.sav_wide(:,p.nyP,iyF,iyb)./xgrid.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',1);
    hold on;
    plot(sgrid.short,ones(p.nx,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([p.borrow_lim p.xmax]);
    title('Savings Policy Function s/x');

    % consumption policy function: zoomed in
    subplot(2,4,3);
    plot(xgrid.orig_wide(:,1,iyF),basemodel.con_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.con_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',2);
    grid;
    xlim([0 4]);
    title('Consumption: Zoomed');

     % savings policy function: zoomed in
    subplot(2,4,4);
    plot(xgrid.orig_wide(:,1,iyF,iyb),basemodel.sav_wide(:,1,iyF,iyb)./xgrid.orig_wide(:,1,iyF,iyb),'b-',xgrid.orig_wide(:,p.nyP,iyF,iyb),basemodel.sav_wide(:,p.nyP,iyF,iyb)./xgrid.orig_wide(:,p.nyP,iyF,iyb),'r-','LineWidth',2);
    hold on;
    plot(sgrid.short,ones(p.nx,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([0 4]);
    title('Savings (s/x): Zoomed');

     % gross income distribution
    subplot(2,4,5);
    b = bar(income.ysort,income.ysortdist);
    b.FaceColor = 'blue';
    b.EdgeColor = 'blue';
    grid;
    xlim([0 10]);
    title('Gross Income PMF');
    
    % create histogram of asset holdings
    binwidth = 0.25;
    bins = 0:binwidth:p.xmax;
    values = zeros(p.xmax+1,1);
    ibin = 1;
    for bin = bins
        if bin < p.xmax
            idx = (basemodel.sav_longgrid_sort>=bin) & (basemodel.sav_longgrid_sort<bin+binwidth);
        else
            idx = (basemodel.sav_longgrid_sort>=bin) & (basemodel.sav_longgrid_sort<=bin+binwidth);
        end
        values(ibin) = sum(basemodel.SSdist_sort(idx));
        ibin = ibin + 1;
    end

     % asset distribution
    subplot(2,4,6);
    b = bar(bins,values);
    b.FaceColor = 'blue';
    b.EdgeColor = 'blue';
    grid;
    xlim([-0.4 10]);
    ylim([0 1]);
    title('Asset Histogram');

     % simulation convergence
    if p.Simulate == 1
        subplot(2,4,7);
        plot(1:p.Tsim,assetmeans,'b','LineWidth',2);
        grid;
        xlim([0 p.Tsim]);
        title('Mean savings (sim)');
    end


end