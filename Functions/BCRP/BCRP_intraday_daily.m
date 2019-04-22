function bcrp_perf = BCRP_intraday_daily(x_daily,x_intraday)            

    tmin0 = 2;
    T = dimdaily(3); %total number of time periods in the data

    if T < tmin0
        error('pattern:offline','Not enough Data L*K>T');
    end

    % compute returns for periods 1 to tmin0 (before all trading
    % starts)
    % Extract closing prices to compute returns
    Pc_daily = reshape(x_daily(:,1,1:tmin0-1),size(x_daily(:,1,1:tmin0-1),1),size(x_daily(:,1,1:tmin0-1),3)); %(stocks x time)

    % Method 1: Compute arithmetic returns for closing s and returns
    %FIX ME!!
    %ret = [zeros(1,size(Pc,1)+1);[diff(Pc')./Pc(:,1:end-1)',diff(STEFI(1:t,1))./STEFI(1:(t-1),1)]]; %add row of zeros at the start of vector to account for 1st period ret

    % Method 2: Compute log returns
    %FIX ME!!
    %ret = [zeros(1,size(Pc,1)+1);[diff(log(Pc),1,2)',diff(log(STEFI(1:t,1)))]]; %add row of zeros at the start of vector to account for 1st period ret
    %ret = ret(:,[liquidstocks(:,1:noofportfoliostocks) end]);

    % Method 3: Compute price relatives
    ret_all_daily = [ones(1,size(Pc_daily,1)+1);[transpose(Pc_daily(:,2:end)./Pc_daily(:,1:end-1)),STEFI(2:tmin0-1,1)./STEFI(1:tmin0-1-1,1)]]; %only compute returns up to befre tmin0 as tmin0

    % set controls for initial time period (t=1) as no trading takes place throughout this day
    h(:,:,1:uniquedaytr(2)) = zeros;%(size(h(:,:,1)));
    b(1:uniquedaytr(2),:) = zeros;

    % implement the sequential learning for next time t
    for t = tmin0:T   %day loop
       %initialise intraday returns for t-th day 
       ret_all_intra = ones(1,size(Pc_daily,1)+1);

       STEFI_intra = repmat(STEFI(t,1),uniquedaytr(t+1) - uniquedaytr(t),1);

       for t_intra = (uniquedaytr(t)+1):uniquedaytr(t+1)-1 %intraday loop: +2 and not -1 from uniquedaytr(t+1) since we want the first entry of each day in h, S, etc to be the daily trade
           Pc_intra = reshape(x_intraday(:,1,t_intra-1:t_intra),size(x_intraday(:,1,t_intra-1:t_intra),1),size(x_intraday(:,1,t_intra-1:t_intra),3)); %closing prices time t-1 to t  
           ret_all_intra = [ret_all_intra; [transpose(Pc_intra(:,end)./Pc_intra(:,end-1)),STEFI_intra(t_intra-uniquedaytr(t)+1,1)./STEFI_intra(t_intra-uniquedaytr(t),1)]];
       end
       %Daily learning 
       Pc_daily = reshape(x_daily(:,1,t-1:t),size(x_daily(:,1,t-1:t),1),size(x_daily(:,1,t-1:t),3)); %closing prices time t-1 to t                
       ret_all_daily = [ret_all_daily; [transpose(Pc_daily(:,end)./Pc_daily(:,end-1)),STEFI(t,1)./STEFI(t-1,1)]];
    end   
            
    % now that we have the returns, construct random portfolios and assess
    % performance to get max hindsight portfolio
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


end