function x = constraint_calibrator(beta_in, p)
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	p.beta0 = beta_in;
	results = main(p);

	x = results.direct.wealth_lt_1000 - 0.23;

	fprintf("\n\n---- For beta = %f, P(assets<$1000) = %f ----\n\n",...
		p.beta0, results.direct.wealth_lt_1000)

end