function out = RSIstrat(Data, n1)
    %% Masters: RSI Trading Strategy function
    %
    % Function to calculate the relative strength index of a data set and
    % implement the RSI rule to generate buy/sell/hold signals.
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 is the look-back period for the indicator
    
    % Redefine the function inputs
    n1 = Data{2};
    Data = Data{1};
    
    % Extract the required data- this extracts only the closing prices for
    % periods 1:t as required by 
    Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    
    signal = zeros(size(Datanew,1),1);
    RSI = zeros(size(Datanew,2),size(Datanew,1));
    % Calculate RSI for current and previous period for stock i
    for i = 1:size(Datanew,1) %loop over stocks
        RSI(:,i) = rsindex(Datanew(i,:)',n1);
        %RSIt1 = rsindex(Datanew(i,end-n1:end)',n1); 
        RSIt1 = RSI(end,i);
        RSIt0 = RSI(end-1,i);
        % Compute signals for rule
        if RSIt0 <= 30 && RSIt1 > 30
         signal(i,1) = 1;
        elseif RSIt0 >= 70 && RSIt1 < 70
         signal(i,1) = -1;
        else
         signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
    out{2} = RSI;
end
