function out = MARSI(Data, n1)
    %% Masters: RSI Trading Strategy function
    %
    % Function to calculate the relative strength index of a data set and
    % implement the RSI rule to generate buy/sell/hold signals.
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 is the look-back period for the indicator
    %https://pdfs.semanticscholar.org/f2c8/ee4e0b6f81eeba97413f5a2b4b120258618e.pdf
    %MARSI stands for Moving Average Relative Strength Index.
%     MARSI is an indicator used to smooth out the action of RSI
%     indicator. This indicator is the same a RSI except that we
%     calculate a Simple Moving Average (SMA) for the calculated
%     RSI indicator. Instead of buying or selling when RSI crosses
%     the thresholds, Buy and sell signals are generated when the
%     moving average crosses above or below the threshold levels.
% 
%     The parameters of the MARSI indicator are: The number of
%     days for look-back period, the lower and upper thresholds and
%     the number of days to average the RSI.

    
    % Redefine the function inputs
    n1 = Data{2};
    Data = Data{1};
    
    % Extract the required data- this extracts only the closing prices for
    % periods 1:t as required by 
    Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    
    RSI = nan(size(Datanew,2),size(Datanew,1));
    MARSI = RSI;
    
    signal = zeros(size(Datanew,1),1);
    % Calculate RSI for current and previous period for stock i
    for i = 1:size(Datanew,1) %loop over stocks
        RSI(:,i) = rsindex(Datanew(i,:)',n1)';
        
        %SMA of RSI for last n1 periods
        MARSI(:,i) = SMA(RSI(:,i),n1);
        
        % current and previous period RSI values
        MARSIt1 = MARSI(end,i);
        MARSIt0 = MARSI(end-1,i);
        
        % Compute signals for rule
        if MARSIt0 <= 30 && MARSIt1 > 30
         signal(i,1) = 1;
        elseif MARSIt0 >= 70 && MARSIt1 < 70
         signal(i,1) = -1;
        else
         signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
    out{2} = MARSI;
end