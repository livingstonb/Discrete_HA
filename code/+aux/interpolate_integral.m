function interpolant = interpolate_integral(gridValues, integrandValues, pmf, is_sorted)
	% Returns an interpolant that interpolates to find the value of
	% int_0^{epsilon} values(a)g(a)da for a given epsilon.
	%
	% Evaluate the above integral by calling interpolant(epsilon).
	%
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	if nargin < 4
		is_sorted = false;
	end

	if ~is_sorted
		sortedInputs = sortrows([gridValues(:) integrandValues(:) pmf(:)]);
		gridSorted = sortedInputs(:,1);
		integrandSorted = sortedInputs(:,2);
		pmfSorted = sortedInputs(:,3);
	else
		gridSorted = gridValues;
		integrandSorted = integrandValues;
		pmfSorted = pmf;
	end

	integralValues= cumsum(integrandSorted .* pmfSorted);

	dsupport = pmfSorted > 1e-7;
	integralValues = integralValues(dsupport);
	gridSorted = gridSorted(dsupport);

	[gridUnique, uniqueInds] = unique(gridSorted,'last');
	integralUnique = integralValues(uniqueInds);

	interpolant = griddedInterpolant(gridUnique, integralUnique,'linear');

	xmin = gridUnique(1);
	xmax = gridUnique(end);
	int0 = integralUnique(1);
	int1 = integralUnique(2);
	interpolant = @(x) adjust_interpolant(interpolant, x, xmin, xmax, int0, int1);
end

function vals_adj = adjust_interpolant(interpolant0, x, xmin, xmax, int0, int1)
	vals_adj = interpolant0(x);
	vals_adj(x<=xmin) = int0;
	vals_adj(x>=xmax) = int1;
end