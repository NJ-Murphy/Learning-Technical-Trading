function [stats,stats_pnl] = expert_perf_stats(Sh,PnL_experts,params,strategy_names)
    % Sh = xlsread("C:/Users/nicjm/Documents/MATLAB/Masters/Experts_perf/stocks20_strats16_02-Jan-2018-27-Mar-2018_60_PnLexperts.xlsx"); 
    % params = xlsread("C:/Users/nicjm/Documents/MATLAB/Masters/Experts_perf/stocks20_strats16_02-Jan-2018-27-Mar-2018_60_parameters.xlsx"); 

    %% Analyse stats of terminal wealths (S)
    %rank the experts terminal wealths
    [~,rank] = sort(Sh(end,:));

    uc1 = unique(params(:,2)) ;
    mc2 = accumarray(params(:,2), rank, [], @mean) ;
    %newData = [uc1, mc2(uc1)];

    [~,~,c] = unique(params(:,2));

    stats = zeros(5,length(uc1));
    stats(1,:) = accumarray(c,rank,[],@mean); %[a, accumarray(c,rank,[],@mean)];
    stats(2,:) = accumarray(c,Sh(end,:),[],@mean);
    stats(3,:) = accumarray(c,Sh(end,:),[],@std);
    stats(4,:) = accumarray(c,Sh(end,:),[],@min);
    stats(5,:) = accumarray(c,Sh(end,:),[],@max);

    %% Analyse stats of PnL's
    stats_pnl = zeros(length(strategy_names),4);
    
    for j = 1:length(strategy_names)
        currentstrat = PnL_experts(:,c==j);
        stats_pnl(j,1) = mean(currentstrat(:));
        stats_pnl(j,2) = std(currentstrat(:));
        stats_pnl(j,3) = min(currentstrat(:));
        stats_pnl(j,4) = max(currentstrat(:));
    end
    %% Best and worst expert in terms of terminal S and which strategy
   % bestexp = strats(params(rank==1,2));
   % worstexp = strats(params(rank==size(params,1),2));
end