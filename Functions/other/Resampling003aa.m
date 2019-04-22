function timebars = Resampling003aa(frequency,tradedata)
%% Resampling data function for Tick Data
%
% Author: N.J. Murphy
%
% Uses: Follows from Resampling002 and is edited to work with Bloomberg
% tick data dowmloaded from UCT Bloomberg terminal. 

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Process the Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % Step 1: find indices for unique days, start of days and end of days
  [yr,mt,dy,~,~,~] = datevec(tradedata(:,1));
  tradedate = datenum(yr,mt,dy);
  [~,uniquedaytr] = unique(tradedate,'first');  %find the indices for each day 
  
  % Step 2: Create an empty cell array to store the aggregated trades into 1 minute time bars
  % and find the date vectors of all the trades. Generate the time bars for
  % a single day from 9:15 until 16:30 which is when closing auctions take place. 
  formatIn = 'HH:MM:SS';  %format of time needed for various functions
  [yrs,mts,dys,hrs,mns,scs] = datevec(tradedata(:,1)); %find year, month, day, hour, minute and second of each trade
  [tradedatetimes] = datenum(yrs,mts,dys,hrs,mns,scs);  %input trade date-times into vector
  minsinday = 50 + 6*60 + 35; %no. minutes per day in continuous trading (excl.first 10min and last 20min) from 9:15am to 16:30pm GMT
  nooftimebars = minsinday/frequency;  %number of time bars per day
  timebars = cell(nooftimebars*length(uniquedaytr),9); %create array to store time bars and corresponding data
  %times = cellstr(datestr(9/24:1/1440:17/24,formatIn));  %create a character string of 1 minute time bars- fill in beginning of time-bar
  times = cellstr(datestr(9/24+(1/6)/24:frequency/1440:(16+3/6)/24,formatIn));  %create a character string of 'frequency'-size minute time bars- fill in beginning of time-bar
 
  % Step 3: Fill in the dates and times of each time bar over the entire
  % period into the empty cell array
  ind = 0;  %create a temporary index variable
  for i = 1:length(uniquedaytr) %loop over days
    for j = 1:nooftimebars  %loop from beginning to end of day
      timebars{j+(i-1)*nooftimebars,1} = tradedata{j+uniquedaytr(i),1}; %fill in dates of time bar
      timebars(j+(i-1)*nooftimebars,2) = times(j,1);  %fill in times of timeimebars bar
    end
    ind = ind + nooftimebars + 1;  %use index to call times for each day
  end
  [yrbar,mtbar,dybar,~,~,~] = datevec(string(timebars(1:nooftimebars*length(uniquedaytr),1))); %find year, month and day for dates in time bars
  [~,~,~,hrbar,mnbar,scbar] = datevec(string(timebars(1:nooftimebars*length(uniquedaytr),2)));  %find hour, minute and second for dates in time bars
  [tradedatetimebar] = datenum(yrbar,mtbar,dybar,hrbar,mnbar,scbar); %find serial number for dates and times in time-bars

  
  % Step 4: Find the indices of the trades which lie in each of the time
  % bars and sort the trades into the corresponding time bars
  for i = 1:length(uniquedaytr) %loop over days
     for j = 1:nooftimebars  %loop from beginning to end of day
      if j == 1
        [idx,~] = find(tradedatetimes(:,1) <= tradedatetimebar(j+(i-1)*nooftimebars));  %find indices of trades that occur in current time bar
        timebars{j+(i-1)*nooftimebars,3} = tradedata(idx(1:end),:);  %capture all the trades which lie in the jth time bar  
      else
        [idx,~] = find(tradedatetimes(:,1) > tradedatetimebar(j-1+(i-1)*nooftimebars) & tradedatetimes(:,1) < tradedatetimebar(j+(i-1)*nooftimebars));  %find indices of trades that occur in current time bar
        timebars{j+(i-1)*nooftimebars,3} = tradedata(idx(1:end),:);  %capture all the trades which lie in the jth time bar bar
      end
    end
  end

  clearvars dayendtr daystarttr yrbar dybar mtbar scbar hrbar...
      scs mns hrs mts yrs mnbar tradedatetimes times
  
  %Remove first bar of each day as it is not required
  %tradedatetimebar(strcmp(timebars(:, 2),'09:10:00'), :) = [];
  timebars(:,1) = cellstr(datestr(tradedatetimebar,'yyyy-mm-dd HH:MM:SS'));
  timebars(strcmp(timebars(:, 2),'09:10:00'), :) = [];
  
  % Step 5: Compute volume, VWAP, open, high, low, close for each bar
  for i = 1:length(timebars)
    if isempty(timebars{i,3}) == 1
        timebars(i,4:9) = {NaN};
    else
        %stocks other than NPN:
%         timebars{i,4} = (sum(str2double(timebars{i,3}(:,3)).*str2double(timebars{i,3}(:,2))))/sum(str2double(timebars{i,3}(:,3))); %vwap
%         timebars{i,5} = sum(str2double(timebars{i,3}(:,3))); %volume
%         timebars{i,6} = str2double(timebars{i,3}{1,2}); %open
%         timebars{i,7} = max(str2double(timebars{i,3}(:,2))); %high
%         timebars{i,8} = min(str2double(timebars{i,3}(:,2))); %low
%         timebars{i,9} = str2double(timebars{i,3}{end,2});%close
        %NPN:
        timebars{i,4} = (sum(cell2mat(timebars{i,3}(:,3)).*cell2mat(timebars{i,3}(:,2))))/sum(cell2mat(timebars{i,3}(:,3))); %vwap
        timebars{i,5} = sum(cell2mat(timebars{i,3}(:,3))); %volume
        timebars{i,6} = timebars{i,3}{1,2}; %open
        timebars{i,7} = max(cell2mat(timebars{i,3}(:,2))); %high
        timebars{i,8} = min(cell2mat(timebars{i,3}(:,2))); %low
        timebars{i,9} = timebars{i,3}{end,2};%close
    end
  end

  %Remove irrelevant transactions
  timebars(:,2:3) = [];
     
%   % Extract relevant variables: Open, Close, High, Low, Date-Time
%   DateTime = timebars(:,1);
%   P_vwap = cell2mat(timebars(:,2));
%   Vol = cell2mat(timebars(:,3));
%   P_o = cell2mat(timebars(:,4));
%   P_c = cell2mat(timebars(:,7));
%   P_h = cell2mat(timebars(:,5));
%   P_l = cell2mat(timebars(:,6));
end