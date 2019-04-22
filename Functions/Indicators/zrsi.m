function out = rsi(data,period)
% Function to calculate the Relative Strength Index of a data set
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
% 'period' is the length of the Wilder Smoothing window.
%
% Example:
% out = rsi(data,period)
%

% Error check
if nargin ~=2 
    error([mfilename,' requires 2 inputs.']);
end
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The first input to ',mfilename,' must be a vector.']);
end
if numel(period) ~= 1
    error('The Wilder smoothing period must be a scalar.');
end
if length(data) < period+1
    error('The data set must be at least 1 element longer than the requested RSI period.');
end

% calculate the up and down data changes
dd = diff(data);
uc = dd;
uc(uc<0)=0;
dc = dd;
dc(dc>0)=0;
dc = -dc;
% perform Wilder Smoothing
wuc = wildersmoothing(uc,period);
wdc = wildersmoothing(dc,period);
% calculate the RSI (taking account of nan's in wuc and wdc)
out = [nan*ones(period-1,1); 100-(100./(1+wuc(period-1:end)./wdc(period-1:end)))];