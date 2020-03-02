function mpcs_a = collapse_mpcs(mpcs_states, pmf, pmf_a)
	na = size(mpcs_states, 1);
	mpcs_states = reshape(mpcs_states, na, []);

	mpcs_a = sum(mpcs_states .* pmf, 2)...
		./ pmf_a;

	pmf_a_small = pmf_a < 1e-8;
	mpcs_a(pmf_a_small) = mean(mpcs_states(pmf_a_small,:), 2);
end