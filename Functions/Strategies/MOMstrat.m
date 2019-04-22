function out = MOMstrat(Data, n1)
    %% Masters: Momentum Trading Strategy function
    %
    % Function to calculate the momentum of a data set
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
    
    % Resize the input data as only closing prices are needed 
    Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    clearvars Data
    

    %% Compute MOM and EMA(MOM)

    % Initialise variables
    mom = nan(size(Datanew,1),size(Datanew,2));
    %noofobserv = size(Datanew,2);
    signal = nan(size(Datanew,1),1);
    
    % Loop over stocks
    for i = 1:size(Datanew,1)
        % EMA of Long Period
        mom(i,:) = tsmom(Datanew(i,:)',n1);
        [EMAMOM, status] = ema(mom,n1,length(mom));
        if ~status
           signal = zeros(size(signal));
           %out{1} = signal;
           break
        end
        
        % Store the t-1 (previous) and t (current) values of MACD and MACDS
        momt0 = mom(i,end-1);
        momt1 = mom(i,end);
        emamomt1 = EMAMOM(end);
        
        % Generate the buy/sell/hold signals
        if momt0 <= emamomt1 && momt1 > emamomt1
            signal(i,1) = 1;
        elseif momt0 >= emamomt1 && momt1 < emamomt1
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end        
    
    out{1} = signal;
end