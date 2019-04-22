function SMAout = SMA(data,n)
% Function to calculate the simple moving average of a data set
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
% 'n' is the number of periods over which to calculate the average
%
% Example:
% SMAout = sma(data,n)
%

% Error check
if nargin ~= 2
    error([mfilename,' requires 2 inputs.']);
end
[r,c]=size(data);
if ~(r==1 || c==1)
    error(['The data input to ',mfilename,' must be a vector.']);
end
if (numel(n) ~= 1) || (mod(n,1)~=0)
    error('The n must be a scalar integer.');
end
if length(data) < n
    error('The length of the data must be at least the specified ''n''.');
end

% calculate the SMA
SMAout = filter(ones(1,n),n,data);
SMAout(1:n-1) = nan; % these are just the filter buffer filling

end