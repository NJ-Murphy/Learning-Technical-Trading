function [ret,b] = BCRP_daily(x,n_portfolios)
  
    n_assets = size(x,2);

    % generate a set of random porttfolios
    portfolios_weights = gen_port(n_portfolios, n_assets);

     % calculate terminal wealth of sample portfolios
     portfolios_wealth = zeros(size(portfolios_weights,1),1);
     period_wealth = zeros(size(portfolios_wealth,1), size(x,1));

     for i = 1:size(portfolios_weights,1)
        period_wealth(i,:) = x * portfolios_weights(i,:)';
        port_ret = cumprod(period_wealth(i,:));
        portfolios_wealth(i) = port_ret(end);
     end
  
    % Find max CRP portfolio wealth at time T
    [~,BCRP_max] = max(portfolios_wealth);
    % BCRP weights of best performng portfolio
    b = portfolios_weights(BCRP_max,:)';

    % BCRP performance
    S = cumprod(x*b);
    ret = [1; S];

end

