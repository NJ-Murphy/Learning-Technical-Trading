function out = WilliamsR(Data, n1)
%% Masters: Williams Percent Range Trading Strategy function
    %
    % Function to calculate the fast stochastic trading rule for a data set
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 is the number of periods over which to calculate the simple moving
    % averages.
    % n is the period of the indicator moving average
    
    % Redefine the function inputs
    n1 = Data{2};
    Data = Data{1};
    
    % Extract the required data- this extracts close, high and low
    %Datanew = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    hi = Data(:,3,:);
    lo = Data(:,4,:);
    cl = Data(:,1,:);
    
    hi = reshape(hi(:,1,:),size(hi(:,1,:),1),size(hi(:,1,:),3))';
    lo = reshape(lo(:,1,:),size(lo(:,1,:),1),size(lo(:,1,:),3))';
    cl = reshape(cl(:,1,:),size(cl(:,1,:),1),size(cl(:,1,:),3))'; 
   
    %% Williams %R
    
        % call willpctr
        wpcr = willpctr(hi, lo, cl, n1);
        willt1 = wpcr(end-1);
        willt2 = wpcr(end); 
        
        % Initialise variables
       signal = nan(size(cl,2),1);

    for i = 1:size(cl,2)     

        % Generate the buy/sell/hold signals
        if willt1  >=  -20  && willt2  <  -80
            signal(i,1) = 1;
        elseif willt1  <=  -20  && willt2  >  -80
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
    
end