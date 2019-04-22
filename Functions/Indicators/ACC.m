function out = ACC(Data, n1)
%% Masters: Acceleration Trading Strategy function
    %
    % Function to calculate the acceleration of a data set
    % 'Data' is the matrix to operate on.  The first element is assumed to be
    % the oldest data.
    %
    % n1 is the number of periods over which to calculate the price change
    
    % Redefine the function inputs
    %n1 = Data{2};
    %Data = Data{1};
    
    % Extract the required data- this extracts only the closing prices for
    % periods 1:t as required by 
    %Pc = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3));
    
    mom = MOM(Data,n1);

    ACC = mom(end) - mom(end-n1);
    out = ACC;
end