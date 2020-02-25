
clear
close all
load('/home/brian/Documents/temp/chi_tests.mat')

figpath = '/home/brian/Documents/temp/figs';

failures = (chi1 == 0);

chi1 = chi1(~failures);
chi2 = chi2(~failures);
htm = htm(~failures);
fraction_whtm = fraction_whtm(~failures);
chivar = chivar(~failures);

sorted_mat = sortrows([chi1 chi2 htm fraction_whtm chivar], [1 2]);
chi1 = sorted_mat(:,1);
chi2 = sorted_mat(:,2);
htm = sorted_mat(:,3);
fraction_whtm = sorted_mat(:,4);
chivar = sorted_mat(:,5);



chi1vals = sortrows(unique(chi1));
chi2vals = sortrows(unique(chi2));

nchi1 = numel(chi1vals);
nchi2 = numel(chi2vals);

% Plot fraction total HtM vs chi2
fig = figure();
hold on;
legend();
for ichi1 = 1:nchi1
    chi1_val = chi1vals(ichi1);
    idx = (chi1 == chi1_val);
    plot(chi2(idx), htm(idx));
    
    ax = gca;
    ax.Legend.String{end} = sprintf('chi1 = %g', chi1_val);
end
title('Total HtM vs chi2')
xlabel('chi2')
ylabel('P(liquid w <= quart inc / 6)')

saveloc = fullfile(figpath, 'total_htm_vs_chi.jpg');
saveas(gcf, saveloc);

% Plot WHtM / THtM vs chi2
fig = figure();
hold on;
legend();
for ichi1 = 1:nchi1
    chi1_val = chi1vals(ichi1);
    idx = (chi1 == chi1_val);
    plot(chi2(idx), fraction_whtm(idx));
    
    ax = gca;
    ax.Legend.String{end} = sprintf('chi1 = %g', chi1_val);
end
title('WHtM/HtM vs chi2')
xlabel('chi2')
ylabel('P(Wealthy HtM) / P(Total HtM)')

saveloc = fullfile(figpath, 'fraction_whtm_vs_chi.jpg');
saveas(gcf, saveloc);


% Plot chivar vs chi2
fig = figure();
hold on;
legend();
for ichi1 = 1:nchi1
    chi1_val = chi1vals(ichi1);
    idx = (chi1 == chi1_val);
    plot(chi2(idx), chivar(idx));
    
    ax = gca;
    ax.Legend.String{end} = sprintf('chi1 = %g', chi1_val);
end
title('adjust cost multiplier vs chi2')
xlabel('chi2')
ylabel('$\chi_1 ^{-\chi_2} / (1+\chi_2)$', 'interpreter', 'latex')

saveloc = fullfile(figpath, 'chivar.jpg');
saveas(gcf, saveloc);
