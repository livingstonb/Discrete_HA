clear
close all

load('/home/brian/Documents/income.mat')
ydist = income.ymatdist(:,2);
ydist = ydist(:) / sum(ydist(:));
ygrid = income.ymat(:,2);
[bins, vals] = aux.smoothed_histogram(ygrid, ydist, 20, 0.2);
% histogram('BinEdges', bins, 'BinCounts', vals)

%% Percentiles of income dist
iyP = 2;

yvals = reshape(income.ymat(iyP,:), [], 1);
yprobs = reshape(income.ymatdist(iyP,:), [], 1);
yprobs = yprobs / sum(yprobs(:));

sorted_input = sortrows([income.ymat(:) income.ymatdist(:)]);
ymat = sorted_input(:,1);
ymatdist = sorted_input(:,2);

% Support of the distribution
support = ymatdist > 1e-9;
ymat_support = ymat(support);
pmf_support = ymatdist(support);

% cdf(y) over the support of pmf(y)
tmp = cumsum(ymatdist);
tmp = tmp(support);
cdf_support = tmp / tmp(end);

plot(ymat_support, cdf_support)