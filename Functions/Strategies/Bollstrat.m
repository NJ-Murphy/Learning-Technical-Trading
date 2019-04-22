function out = Bollstrat(Data, n1)
    %% Masters: Bollinger Trading Strategy function
    %
    % Function to calculate the upper, lower and middle bollinger bands of a
    % data set and implement the Bollinger rule to generate buy/sell/hold signals.
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 is the look-back period for the indicator
    
    % Redefine the function inputs
    n1 = Data{2};
    Data = Data{1};
    
    % Extract the required data- this extracts only the closing prices for
    % periods 1:t as required by Bollinger Bands strategy
    Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    
    % Initialise variables
    uppr = nan(size(Datanew,1),size(Datanew,2))';
    mid = nan(size(Datanew,1),size(Datanew,2))';
    lowr = nan(size(Datanew,1),size(Datanew,2))';
    signal = nan(size(Datanew,1),1);
    
    % Calculate RSI for current and previous period for stock i
    for i = 1:size(Datanew,1) %loop over stocks
        [mid(:,i), uppr(:,i), lowr(:,i)] = bollinger(Datanew(i,:), n1);
        Bollupt1 = uppr(end,i);
        Bolldt1 = lowr(end,i);
        Pt0 = Datanew(i,end-1);
        Pt1 = Datanew(i,end);
        
        % Compute signals for rule
        if Bolldt1 <= Pt0 && Pt1 >=Bollupt1
         signal(i,1) = 1;
        elseif Bolldt1 >= Pt0 && Pt1 > Bollupt1
         signal(i,1) = -1;
        else
         signal(i,1) = 0;
        end
    end
 
    out{1} = signal;
    out{2} = lowr;
    out{3} = mid;
    out{4} = uppr;
    
end
