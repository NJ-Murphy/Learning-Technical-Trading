function out = Faststrat(Data, n1)
%% Masters: Fast Stochastic Trading Strategy function
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
    
    %this can be made faster by replacing below in above step
    hi = reshape(hi(:,1,:),size(hi(:,1,:),1),size(hi(:,1,:),3))';
    lo = reshape(lo(:,1,:),size(lo(:,1,:),1),size(lo(:,1,:),3))';
    cl = reshape(cl(:,1,:),size(cl(:,1,:),1),size(cl(:,1,:),3))'; 
     
    %% Fast%K and Fast%D
    
    % call Fast%K and Fast%D function
    [pctk, pctd] = fpctkd(hi, lo, cl, n1, 3);
    FastK1 = pctk(end-1);
    FastD2 = pctd(end); 
    %FastD1 = pctd(end-1);
    FastK2 = pctk(end); 
        
    % Initialise variables
    signal = nan(size(hi,2),1);
    
    for i = 1:size(cl,2)           
        % Generate the buy/sell/hold signals
        if FastK1  <=  FastD2  && FastK2 > FastD2
            signal(i,1) = 1;
        elseif FastK1  >=  FastD2  && FastK2 < FastD2
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
    
end