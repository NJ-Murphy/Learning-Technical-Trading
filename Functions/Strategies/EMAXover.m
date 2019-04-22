function out = EMAXover(Data, n1, n2)

    %% Masters: Exponential Moving Average Crossovers function
    %
    % Function to calculate the Exponential Moving Average Crossovers of a data set
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 and n2 are the number of periods over which to calculate the exp moving
    % averages that are subtracted from each other.
    % n is the period of the indicator moving average
    %https://pdfs.semanticscholar.org/f2c8/ee4e0b6f81eeba97413f5a2b4b120258618e.pdf
%      EMA is performed by using double
%     EMA crossovers. A buy signal is generated when the shorter
%     moving average (short term moving average) is crossed above
%     the longer moving average (long term moving average)
%     because this represents the beginning of an uptrend.
% 
%     A sell signal is generated when the shorter moving average is
%     crossed below the longer moving average because this
%     represents the beginning of a down trend.
% 
%     One of the famous combinations of moving averages is 20-50
%     crossovers.

    % Redefine the function inputs
    n1 = Data{2};
    n2 = Data{3};
    Data = Data{1};
        
    % Resize the input data as only closing prices are needed 
    Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    clearvars Data
    

    %% Compute MACD and MACD signal for all stocks

    % Initialise variables
    noofobserv = size(Datanew,2);
    signal = nan(size(Datanew,1),1);
    
    % Loop over stocks
    for i = 1:size(Datanew,1)
        
        % EMA of Long Period
        [EMAn2, status] = ema(Datanew(i,:),n2,noofobserv);
        if ~status
            signal = zeros(size(signal));
            out{1} = signal;
            break
            %return
        end
        
        % EMA of Short Period
        EMAn1 = ema(Datanew(i,:),n1,noofobserv);
        
        % current and previous EMA
        eman1t0 = EMAn1(end-1);
        eman2t0 = EMAn2(end-1);
        eman1t1 = EMAn1(end);
        eman2t1 = EMAn2(end);
        
        % Generate the buy/sell/hold signals
        if eman1t1 >= eman2t1 && eman1t0 < eman2t0
            signal(i,1) = 1;
        elseif eman1t1 <= eman2t1 && eman1t0 > eman2t0
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end    
    if any(isnan(signal))
        signal = zeros(size(signal));
    end
    
    out{1} = signal;
    %out{2} = EMAn1;
    %out{3} = EMAn2;
end