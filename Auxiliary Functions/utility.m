function u = utility(risk_aver,c)
	u = zeros(size(risk_aver));
	Ilog_u = (risk_aver == 1);

	u(Ilog_u) = log(c(Ilog_u));
	u(~Ilog_u) = (c(~Ilog_u) .^(1-risk_aver(~Ilog_u) )-1)./(1-risk_aver(~Ilog_u) );
end