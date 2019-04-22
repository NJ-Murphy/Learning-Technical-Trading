function out = SAR(Data, n1)
%% Edited by N. Murphy
% SAR Parabolic SAR (Stop And Reverse) filter
%
% [SAR,TRND] = SAR(PRC) The Parabolic SAR for price timeseries PRC. PRC can 
%   be a NxM matrix of returns for M stocks and N periods. TRND is the 
%   trend time-series as +1 for up-trend and -1 for down-trend.
%
% [SAR,TRND] = SAR(HI,LO,CL,OP) for high (HI), low (LO), close (CL) and open (OP)
%   Here HI, LO, CL and OP are NxM matrices for M stocks and N periods.
%
% [SAR,TRND] = SAR(HI,LO,CL,OP,ACC) Acceleration increment ACC is by default 0.02.
%
% [SAR, TRND] = SAR(HI,LO,CL,OP,ACC,EP) Extreme-Point EP initialised with price.
%
% This is an indicator developed by J. Welles Wilders Jr. with a trailing stop. SAR can
% never go backwards, as a trend continues the SAR converges on the price
% until is crosses the price.
%
% SAR(i) = SAR(i-1) + ACC [EPRICE(i-1)-SAR(i-1)]
%
% 1. SAR(i-1) is value of indicator
% 2. ACC acceleration factor initialised at 0.02
% 3. EPRICE(i-1) is the highest(lowest) price for the observation period

% TBDL: include the high and low price data for the modification rules

% $Id: sar.m 82 2011-08-10 09:04:07Z tgebbie $

% initial values
% Redefine the function inputs
%n1 = Data{2};
Data = Data{1};

%op = Data(:,2,:);
hi = Data(:,3,:);
lo = Data(:,4,:);
cl = Data(:,1,:);

% reshape input vectors
% op = reshape(op(:,1,:),size(op(:,1,:),1),size(op(:,1,:),3))';
hi = reshape(hi(:,1,:),size(hi(:,1,:),1),size(hi(:,1,:),3))';
lo = reshape(lo(:,1,:),size(lo(:,1,:),1),size(lo(:,1,:),3))';
cl = reshape(cl(:,1,:),size(cl(:,1,:),1),size(cl(:,1,:),3))';

ma = 0.20; % maximum allowed acceleration
a0   = 0.02; % default acceleration increment
%ep = []; % extreme-point is initialised as empty

% initialise the data variables
%sar = nan(size(cl));
% initialise filter
sar = cl;

% initialise the high(low) price
ep = zeros(size(cl,2)); 

% initialise the acceleration
a = a0*ones(1,size(cl,2));

% initialise the lastHigh and lastLow counters.
lastHigh = hi(1,:);
lastLow = lo(1,:);

% initialise the level
trend = ones(size(cl)); % start long

% apply the filter
for j = 1:size(cl,2) %stock loop
    for i = 2:size(cl,1)-1 %time loop
        % update price trends and extreme-points
        if hi(i,j)>lastHigh(j) && trend(i,j)>0
            % postive trend continues (update acceleration)
            a(j) = min(a(j)+a0,ma);
            % update the lastHigh
            lastHigh(j) = hi(i,j);
            % new extreme point (on positive trend)
            ep(j) = lastHigh(j);
        elseif lo(i,j)<lastLow(j) && trend(i,j)<0
            % negative trend continues (update acceleration)
            a(j) = min(a(j)+a0,ma);
            % update the lastLow
            lastLow(j) = lo(i,j);
            % new extreme point (on negative trend)
            ep(j) = lastLow(j);
        end
        % compute the parabolic SAR at time i+1 (tomorrow) from i (today)
        % SAR equation --------------------------------------------------------
        sar(i+1,j) = sar(i,j) + a(j) * (ep(j) - sar(i,j));
        % ---------------------------------------------------------------------
        % SAR modification rules (should use high and low)
        if sar(i,j) > lo(i,j) && trend(i,j)>0
            % trend switching modification from up-trend
            % SAR set to last EP
            sar(i+1,j) = lastHigh(j);
            % A set to initial value
            a(j) = a0;
            % new lastLow
            lastLow(j) = lo(i,j);
            % EP set new trend extreme value
            ep(j) = lastLow(j);
            % change the up-trend to a down-trend
            trend(i+1,j) = -1;
        elseif sar(i,j) < hi(i,j) && trend(i,j)<0
            % trend switching modification from down-trend
            % SAR set to last EP
            sar(i+1,j) = lastLow(j);
            % A set to initial value
            a(j) = a0;
            % new lastHigh
            lastHigh(j) = hi(i,j);
            % EP set new trend extreme value
            ep(j) = lastHigh(j);
            % change from down-trend to an up-trend
            trend(i+1,j) = +1; 
        elseif sar(i+1,j) > min(lo(i,j),lo(i-1,j)) && trend(i,j)>0
            % no trend switching modification
            sar(i+1,j) = min(lo(i,j),lo(i-1,j));
            % If Long using daily data: Never move tomorrows applicable 
            % SAR above yesterdays or today's low. If the calculated SAR 
            % is higher than either of these lows then use the lower low 
            % of these two days as the SAR and use this value for SAR 
            % calculation for the next day. 
            trend(i+1,j) = trend(i,j);
        elseif sar(i+1,j) < max(hi(i,j),hi(i-1,j)) && trend(i,j)<0
            % no trend switching modification
            sar(i+1,j) = max(hi(i,j),hi(i-1,j));
            % If Short using daily data: Never move the SAR below the high 
            % of yesterday or today. If the calculated SAR is lower than 
            % either of these values then use the higher high of these as 
            % the SAR for the day and for the future calculation for the 
            % following day's SAR.
            trend(i+1,j) = trend(i,j);
        else
            % no SAR modification and trend continues
            trend(i+1,j) = trend(i,j);
        end
    end
end

%% Signals
% Initialise variables
signal = nan(size(cl,2),1);

% define buy, sell and hold signals
for i = 1:size(sar,2)
  if (sar(end-1,i)>=cl(end-1,i)) && (sar(end,i)<cl(end,i))
    signal(i,1) = 1;
  elseif (sar(end-1,i)<=cl(end-1,i)) && (sar(end,i)>cl(end,i))
      signal(i,1) = -1;
  else
      signal(i,1) = 0;
  end
end

out{1} = signal;
end