function [vout, status] = ema(vin,lag,observ)

% Number of observations
%observ = size(vin,1);

% Preallocate output
vout   = nan(observ,1);

% Set status
status = 1;

% If the lag is greater than or equal to the number of observations
if lag >= observ
    status = 0;
    return
end

% Calculate the exponential percentage
k = 2/(lag+1);

% Calculate the simple moving average for the first 'exp mov avg' value.
vout(lag) = sum(vin(1:lag))/lag;

% K*vin; 1-k
kvin = vin(lag:observ)*k;
oneK = 1-k;

% First period calculation
vout(lag) = kvin(1)+(vout(lag)*oneK);

% Remaining periods calculation
for i1 = lag+1:observ
    vout(i1) = kvin(i1-lag+1)+(vout(i1-1)*oneK);
end

end
%--------------------------------------------------------------------------