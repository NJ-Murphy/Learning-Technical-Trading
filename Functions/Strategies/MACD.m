function out = MACD(Data, n1, n2)
%, MACDt0)
    %% Masters: MACD Trading Strategy function
    %
    % Function to calculate the moving average convergence/divergence of a data set
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 and n2 are the number of periods over which to calculate the exp moving
    % averages that are subtracted from each other.
    % n is the period of the indicator moving average
    
    % Redefine the function inputs
    n1 = Data{2};
    n2 = Data{3};
    Data = Data{1};
    
    % Set a default for the indicator moving average
    %n = n1-3;%3;
    n = 9;
    if n >= n2
        n = n2-n1;
    end
    %Test EMA -> x=(34350 - Datanew(1,10))*(2/11)+ Datanew(1,10);
    
    % Resize the input data as only closing prices are needed 
    Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    clearvars Data
    
    % Calculate the short and long EMA's
%    EMAn1 = tsmovavg(Datanew,'e',n1);
%    EMAn2 = tsmovavg(Datanew,'e',n2);
    %idx = isnan(MACD);
    
    % Compute MACD =  shortEMA - longEMA at time t and time t-1
   % MACDt1 = EMAn1(:,end) - EMAn2(:,end-1);
   % MACDt0 = EMAn1(:,end-1) - EMAn2(:,end-2);
%    MACD = EMAn1 - EMAn2;
    
    
    % Calculate the MACD signal Line
%   MACDS = tsmovavg(MACD,'e',n);

    %% Compute MACD and MACD signal for all stocks

    % Initialise variables
    MACD = nan(size(Datanew,1),size(Datanew,2));
    MACDS = nan(size(MACD));
    noofobserv = size(Datanew,2);
    signal = nan(size(Datanew,1),1);
    
    % Loop over stocks
    for i = 1:size(Datanew,1)
        
        % EMA of Long Period
        [EMAn2, status] = ema(Datanew(i,:),n2,noofobserv);
        if ~status
           signal = zeros(size(signal));
           %out{1} = signal;
           break
        end
        
        % EMA of Short Period
        EMAn1 = ema(Datanew(i,:),n1,noofobserv);
        
        % MACD
        MACD(i,1:noofobserv) = EMAn1'-EMAn2';
        
        % MACD Signal
        [presignal, status] = ema(MACD(i,~isnan(MACD(i,1:noofobserv))),n,noofobserv-n2+1);
        if ~status
            signal = zeros(size(Datanew,1),1);
            %error('hello %d,%d', n1,n2);
        end
        MACDS(i,1:noofobserv) = [nan(n2-1,1);presignal]';
        clearvars presignal
        
        % Store the t-1 (previous) and t (current) values of MACD and MACDS
        MACDt0 = MACD(i,end-1);
        MACDt1 = MACD(i,end);
        MACDSt1 = MACDS(i,end);
        
        % Generate the buy/sell/hold signals
        if MACDt0 <= MACDSt1 && MACDt1 > MACDSt1
            signal(i,1) = 1;
        elseif MACDt0 >= MACDSt1 && MACDt1 < MACDSt1
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end        
    
    out{1} = signal;
    out{2} = MACDS;
    out{3} = MACD;
end