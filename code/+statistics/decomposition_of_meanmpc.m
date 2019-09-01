function decomp = decomposition_of_meanmpc(p,grids,results)
	% this function decomposes the mean mpc into four main terms

	decomp = struct();
	for ia = 1:numel(p.abars)
        decomp(ia).term1 = NaN;
        decomp(ia).term2 = NaN;
        decomp(ia).term3 = NaN;
        decomp(ia).term4 = NaN;
    end

    doDecomposition = (p.nb==1) && (p.EpsteinZin==0) ...
    	&& (p.bequest_weight==0) && (p.temptation==0) && (numel(p.r)==1);

    if ~doDecomposition
    	return
    end

    % representative agent MPC
    m_ra = p.R * (results.direct.beta*p.R)^(-1/p.risk_aver) - 1;

    % mpcs from model with income risk
    mpcs = results.direct.mpcs(5).mpcs_1_t{1,1}; % mpcs
    meanmpc = results.direct.mpcs(5).avg_s_t(1,1);

    % distribution over full state space
    g = results.direct.adist;

    % distribution condensed to asset grid
    g_nr = results.direct.agrid_dist;

    % mpcs from model without income risk
    mpcs_nr  = results.norisk.mpcs1_a_direct{5};
    meanmpc_nr = mpcs_nr(:)' * g_nr(:);

    % interpolate to get the integral of mpc(a) * g(a) between a = 0 and 0.05
    m0g0 = mpcs(:) .* g(:);
    m0g0 = reshape(m0g0,[p.nx_DST p.nyP*p.nyF*p.nb]);
    m0g0 = sum(m0g0,2);
    cum_m0g0 = cumsum(m0g0);
    m0g0interp = griddedInterpolant(grids.a.vec,cum_m0g0,'linear');

    % interpolate to get the integral of mpc_norisk(a) * g_norisk(a)
    m1g0 = mpcs_nr(:) .* g_nr(:);
    m1g0 = reshape(m1g0,p.nx_DST,[]);
    m1g0 = sum(m1g0,2);
    cum_m1g0 = cumsum(m1g0);
    m1g0interp = griddedInterpolant(grids.a.vec,cum_m1g0,'linear');

    % get interpolant for cumulative dist of g_a
    g_a = sum(reshape(g,p.nx_DST,[]),2);
    ginterp = griddedInterpolant(grids.a.vec,cumsum(g_a),'linear');

    % get interpolant for cumdist of g_norisk_a
    g_nr_a = sum(reshape(g_nr,p.nx_DST,[]),2);
    g_nr_interp = griddedInterpolant(grids.a.vec,cumsum(g_nr_a),'linear');

    % perform decomposition
    for ia = 1:numel(p.abars)
        decomp(ia).term1 = m_ra;

        if p.abars(ia) == 0
            zidx = grids.a.matrix(:) <= p.abars(ia);
            norisk_zidx = grids.a.vec <= p.abars(ia);

            decomp(ia).term2 = (mpcs(zidx) - m_ra)' * g(zidx);
            decomp(ia).term3 = (mpcs_nr(~norisk_zidx) - m_ra)' * g_nr(~norisk_zidx);
            decomp(ia).term4 = mpcs(~zidx)' * g(~zidx) - mpcs_nr(~norisk_zidx)' * g_nr(~norisk_zidx);
        else
            abar = p.abars(ia);
            decomp(ia).term2 = m0g0interp(abar) - m_ra * ginterp(abar);
            decomp(ia).term3 = meanmpc_nr - m1g0interp(abar) - m_ra * (1-g_nr_interp(abar));
            decomp(ia).term4 = (meanmpc - m0g0interp(abar)) - (meanmpc_nr - m1g0interp(abar));
            
        end
    end


end