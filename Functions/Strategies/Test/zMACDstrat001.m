function signal = MACDstrat(Data,n1,n2,n)
%% Masters: MACD Trading Strategy function
%
% [EMA1,EMA2] = macd(data,n1,n2,n)
% Function to calculate the moving average convergence/divergence of a data set
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
%
% n1 and n2 are the number of periods over which to calculate the moving
% averages that are subtracted from each other.
% n is the period of the indicator moving average
%
% If called with one output then it will be a two column matrix containing
% both calculated series.
% If called with two outputs then the first will contain the macd series
% and the second will contain the indicator series.
%
% Example:
% EMA1 = macd(data,n1,n2,n);
% [EMA1,EMA2] = macd(data,n1,n2,n);

% Error check
if (nargin < 1) || (nargin >4)
    error([mfilename,' requires between 1 and 4 inputs.']);
end
[m,n]=size(Data);
if ~(m==1 || n==1)
    error(['The data input to ',mfilename,' must be a vector.']);
end

% set some defaults
switch nargin
    case 1
        n1 = 26;
        n2 = 12;
        n = 9;
    case 2
        n2 = 12;
        n = 9;
    case 3
        n = 9;
end

if (numel(n1) ~= 1) || (numel(n2) ~= 1) || (numel(n) ~= 1)
    error('The period must be a scalar.');
end

% calculate the MACD
EMA1 = EMA(Data,n2)-EMA(Data,n1);
% Need to be careful with handling NaN's in the second calculation
idx = isnan(EMA1);
EMA2 = [EMA1(idx); EMA(EMA1(~idx),n)];
switch nargout
    case {0,1}
        varargout{1} = [EMA1 EMA2];
    case 2
        varargout{1} = EMA1;
        varargout{2} = EMA2;
    otherwise
        error('Too many outputs have been requested.');
end

end