function out = PROCstrat(Data, n1)
%% Masters: Momentum Trading Strategy function
    %
    % Function to calculate the price rate of change of a data set
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 is the number of periods over which to calculate the price change
    
    % Redefine the function inputs
    n1 = Data{2};
    Data = Data{1};
    
    % Extract the required data- this extracts only the closing prices for
    % periods 1:t as required by 
    Pc = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));

    signal = zeros(size(Pc,1),1);
    for i = 1:size(Pc,1)
        proc = prcroc(Pc(i,:), n1);
        proct1 = proc(end-1);
        proct2 = proc(end); 
     
        % Generate the buy/sell/hold signals
        if proct1  <= 0 && proct2 > 0 
            signal(i,1) = 1;
        elseif proct1 >= 0 && proct2 < 0 
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
end