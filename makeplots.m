function makeplots(p,xgrid,sgrid,basemodel,income,sim_results,assetmeans)
    % This function makes plots based on fields from 'direct_results'
    % and 'simulation_results' and 'income' as well as plotting MPCs

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
    binmax = 10;
    [income_bins,income_values] = create_bins(p,binmax,20,income.ysort,...
                                                         income.ysortdist);
    subplot(2,4,5);
    b = bar(income_bins,income_values);
    % b = bar(income.ysort,income.ysortdist);
    b.FaceColor = 'blue';
    b.EdgeColor = 'blue';
    grid;
    xlim([0 binmax]);
    ylim([0 max(income_values)]);
    title('Gross Income Histogram');

    % plot asset histogram
    [asset_bins,asset_values] = create_bins(p,20,20,...
                            basemodel.sav_longgrid_sort,basemodel.SSdist_sort);
    subplot(2,4,6);
    b = bar(asset_bins,asset_values);
    b.FaceColor = 'blue';
    b.EdgeColor = 'blue';
    grid;
    xlim([-0.4 10]);
    % ylim([0 1]);
    ylim([0 max(asset_values)]);
    title('Asset Histogram');

     % simulation convergence
    if p.Simulate == 1
        subplot(2,4,7);
        plot(1:p.Tsim,assetmeans,'b','LineWidth',2);
        grid;
        xlim([0 p.Tsim]);
        title('Mean savings (sim)');
    end

    %% MPCs
    if p.ComputeDirectMPC == 1
        % STILL NEED TO GET MPCs sorted!
        figure(2);
        [mpc1_bins,mpc1_values] = create_bins(p,1,20,basemodel.a_longgrid_sort,...
                                                        basemodel.SSdist_sort);
        b = bar(mpc1_bins,mpc1_values);
        b.FaceColor = 'blue';
        b.EdgeColor = 'blue';
        grid;
        xlim([0 1]);
        ylim([0 max(income_values)]);
        title('One-Period MPC Histogram');
    end

    %% Histogram helper function
    function [bins,values] = create_bins(p,binmax,nbins,vals_sort,dist_sort)
        bins = linspace(0,binmax,nbins);
        binwidth = bins(2) - bins(1);
        values = zeros(nbins,1);
        ibin = 1;
        for bin = bins
            if ibin == 1
                idx = vals_sort < bin + binwidth;
            elseif ibin < nbins
                idx = (vals_sort>=bin) & (vals_sort<bin+binwidth);
            else
                idx = (vals_sort>=bin);
            end
            values(ibin) = sum(dist_sort(idx));
            ibin = ibin + 1;
        end
    end

end