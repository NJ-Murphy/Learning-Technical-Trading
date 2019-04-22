function out = ACCstrat(Data, n1)
%% Masters: Acceleration Trading Strategy function
    %
    % Function to calculate the momentum of a data set
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
    
    % compute the accelaration for n1
    acc = tsaccel(Pc', n1, 0);
    acct1 = acc(end-1,:);
    acct2 = acc(end,:); 
     
    signal = zeros(size(Pc,1),1);
    for i = 1:size(Pc,1)
        % Generate the buy/sell/hold signals
        if acct1(i) + 1  <= 0 && acct2(i) + 1 > 0 
            signal(i,1) = 1;
        elseif acct1(i) + 1  >= 0 && acct2(i) + 1 < 0 
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
end