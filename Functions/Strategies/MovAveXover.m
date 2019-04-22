function out = MovAveXover(Data, n1, n2)
%% Masters: Moving Average Crossover Trading Strategy function
    %
    % Function to calculate the moving average convergence/divergence of a data set
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 and n2 are the number of periods over which to calculate the simple moving
    % averages.
    % n is the period of the indicator moving average
    
    % Redefine the function inputs
    n1 = Data{2};
    n2 = Data{3};
    Data = Data{1};
    
    % Extract the required data- this extracts only the closing prices for
    % periods 1:t as required by 
    Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    
    %% Moving average crossover
    
    % Initialise variables
    signal = nan(size(Datanew,1),1);
    sman1 = nan(size(Datanew,1),size(Datanew,2));
    sman2 = nan(size(sman1));
    
    % Loop over stocks
    for i = 1:size(Datanew,1)
        % Calculate short and long trem moving averages
        sman1(i,:) = SMA(Datanew(i,:),n1);
        sman2(i,:) = SMA(Datanew(i,:),n2);
        
        sman1t0 = sman1(i,end-1);
        sman2t0 = sman2(i,end-1);
        sman1t1 = sman1(i,end);
        sman2t1 = sman2(i,end);
        
        %cross = [NaN; diff(sign(sma_n1-sma_n2))];  % Compute diff of signal

        %ix = abs(cross)==2;
        %cross(ix) = cross(ix)/2;  % Scale transitions to +/-1
        %tr_ma(k) = {cross}; % Store as cell
        
        % Generate the buy/sell/hold signals
        if sman1t1 >= sman2t1 && sman1t0 < sman2t0
            signal(i,1) = 1;
        elseif sman1t1 <= sman2t1 && sman1t0 > sman2t0
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
    out{2} = sman1;
    out{3} = sman2;

end