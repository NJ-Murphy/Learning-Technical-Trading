function out = wildersmoothing(data,period)
% Function to perform Wilder Smoothing on a data set
% 'data' is the vector of values to be smoothed.  The first element is assumed to be
% the oldest data.
% 'period' is the length of the smoothing window.
%
% Example:
% out = wildersmoothing(data,period)
%

% Error check
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The first input to ',mfilename,' must be a vector.']);
end
if numel(period) ~= 1
    error('The Wilder smoothing period must be a scalar.');
end

% perform the filtering
ld = length(data);
if ld < period
    error('The data vector must be at least as long as the required period of smoothing.');
elseif ld == period
    out = mean(data);
else
    out = nan*ones(size(data));
    out(period:end) = filter(1/period,[1 -(period-1)/period],data(period:end),sum(data(1:(period-1)))/period);
end