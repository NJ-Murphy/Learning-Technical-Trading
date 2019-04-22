function out = Slowstrat(Data, n1)
%% Masters: Slow Stochastic Trading Strategy function
    %
    % Function to calculate the slow stochastic trading rule for a data set
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 and n2 are the number of periods over which to calculate the simple moving
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
     
    %% Slow Stochastic
  
        % call Fast%K and Fast%D function
        [fpctk, fpctd] = fpctkd(hi, lo, cl, n1, n1);
        if sum(fpctd>0) <= 3
            signal = zeros(size(cl,2),1);
            out{1} = signal;
            return
        end
        [spctk, spctd] = spctkd(fpctk,fpctd,3);
        SlowK1 = spctk(end-1);
        SlowD2 = spctd(end); 
        %SlowD1 = spctd(end-1);
        SlowK2 = spctk(end); 
        
       % Initialise variables
    signal = nan(size(cl,2),1);
     for i = 1:size(cl,2)   
        % Generate the buy/sell/hold signals
        if SlowK1  <=  SlowD2  && SlowK2 > SlowD2
            signal(i,1) = 1;
        elseif SlowK1  >=  SlowD2  && SlowK2 < SlowD2
            signal(i,1) = -1;
        else
            signal(i,1) = 0;
        end
    end
    
    out{1} = signal;
    
end