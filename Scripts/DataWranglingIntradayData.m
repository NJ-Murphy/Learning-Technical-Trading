%% Masters: Data wrangling for Tick Data
%
% Author: N.J. Murphy
%
% Person Version: Masters-functions-Resampling
%
% 1 Problem Specification:  This script aims to 
%
% 2 Data Specification: 6 months worth of order book data for 16 stocks listed
%   on the JSE Top 40. The specified period is from 01 January 2013 to 03 June
%   2013. 
%   In order to run and publish this script, 2 stocks were ralabelled in order
%   to loop over a 2 week period and 2 stocks from the test data was
%   used.
%
% 3 Data Description:
%   "index","times","type","value","size","condcode"
%
%
% 4 Configuration Control:
%        userpath/MATLAB/Masters
%        userpath/MATLAB/Masters       
%        userpath/MATLAB/Masters/Scripts    
%        userpath/MATLAB/Masters/Functions 
%        userpath/MATLAB/Masters/Data
%        userpath/MATLAB/Masters/html        
%
% 5 Version Control: No current version control
%
% 6 References: None
%
% 7 Current Situation: 
% 8 Future Situation: no future situation
%
% Uses: processes the raw tick data to get it into desired form for
% LearnIntraday script. Removes unneccessary columns and extracts data for chosen
% continuous trading times. Also creates time bars at desired frequency
 

%% 3. Clear Workspace
% This section prepares the workspace for the implementation of the script.
clc;            % clear command window
format long g;  % formating for output on comand window
format compact; % formating for output on comand window
clear %clear variables

%%%%%%%%%%%%%%%%%%%
%% 4. Path Setup %%
% Set the project paths
[projectPath] = pwd;
% Add the project path so that the script sees the necessary functions
addpath(projectPath);
userpathstr = userpath;
userpathstr = userpathstr(~ismember(userpathstr,';'));
projectpath = 'Masters';   %path for personal computer
addpath(fullfile(userpathstr,projectpath,'Functions'));
addpath(fullfile(userpathstr,projectpath,'Functions','Strategies'));
addpath(fullfile(userpathstr,projectpath,'Functions','Indicators'));
addpath(fullfile(userpathstr,projectpath,'Scripts'));
addpath(fullfile(userpathstr,projectpath,'Data\Bloomberg\tick'));

% Paths for data
datapath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Data\Bloomberg\tick';
datadirec = dir(fullfile(datapath));
datafilenames = {datadirec(~[datadirec.isdir]).name}';
writepath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Data\Bloomberg\tick\Processed2\';

% choose time bar size
barsize = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Load Data and Process the Data %%

%extract data from .csv file
for l = 1:length(datafilenames)
  % process and store the raw data to get desired format
  [~, ~, rawdata] = xlsread(datafilenames{l}); %load csv

  xc=strrep(rawdata(:,1),'"TRADE",','');
  xc=strrep(xc,',"IP"','');
  xc=strrep(xc,',"UC"','');
  xc=strrep(xc,',"AT"','');
  xc=strrep(xc,',"XO"','');
  xc=strrep(xc,',"OP,PD"','');
  xc=strrep(xc,',"OP,L"','');
  xc=strrep(xc,',"LT,L"','');
  xc=strrep(xc,',"LC,L"','');
  xc=strrep(xc,',"LT,PD"','');
  xc=strrep(xc,',"OP"','');
  xc=strrep(xc,',"LT"','');
  xc=strrep(xc,',"AC"','');
  xc=strrep(xc,',"OC"','');
  xc=strrep(xc,',"PX"','');
  xc=strrep(xc,',"OC"','');
  xc=strrep(xc,',"*X"','');
  xc=strrep(xc,',"PF"','');
  xc=strrep(xc,',"BT"','');
  xc=strrep(xc,',"BK"','');
  xc=strrep(xc,',"BK,L"','');
  xc=strrep(xc,',"BK,PD"','');
  xc=strrep(xc,',"PC,L"','');
  xc=strrep(xc,',"CF"','');
  xc=strrep(xc,',"CF,L"','');
  xc=strrep(xc,',"CF,PD"','');
  
  %rawdatacell = rawdata(cellfun(@(s)isempty(regexp(s,'""')),rawdata)); %remove rows which contained ,L (closing auction values)
  rawdatacell = regexp(xc,',','split'); %split comma delimited values into separate columns
   for i = 1:size(rawdatacell,1)
       if i == 622867
           disp('j')
       end
     temp(i,:) = cellstr(rawdatacell{i,:});
   end
  
  appendrawcell = vertcat(rawdatacell{:}); %reappend cell array
  appendrawcell(:,1) = []; %remove index first row, column describing type of data (trade) 
%   for i = 1:size(rawdatacell,1)
%       if size(rawdatacell{i,1},1) ~= 1 && size(rawdatacell{i,1},2)~= 6
%           disp('j')
%       end
%   end

% convert times to SAST from GMT: add 2 hours to each entry
  %%% --> this works perfectly if you want the result to be a cell array:
  appendrawcell(:,1) = cellstr(datestr(datetime(appendrawcell(:,1))+hours(2),'yyyy-mm-dd HH:MM:SS'));
  tradedata = appendrawcell;
  %%% --> if numeric output needed then use this: converting to numeric
  %%% increases speed ten fold
%   tradedata = zeros(size(appendrawcell));
%   tradedata(:,1) = datenum(datetime(appendrawcell(:,1))+hours(2));
%   for i = 1:size(tradedata,1)
%    tradedata(i,2) = str2num(cell2mat(appendrawcell(i,2)));
%    tradedata(i,3) = str2num(cell2mat(appendrawcell(i,3)));
%   end
  
  % Get the date-time of the trades
  [~,~,~,hr,mn,sc] = datevec(tradedata(:,1)); %find hour, minute and second
  %[tradedatetime] = datenum(yr,mt,dy,hr,mn,sc);  %create vector of dates and times 
  
  % Remove all events between before 09h10 and after 16h50
  removalindextime(1:size(((hr<9) | (hr==9&mn<=10) | (hr==16&mn>=35) | (hr==17)),1)) = ((hr<9) | (hr==9&mn<=10) | (hr==16&mn>=35) | (hr==17)); %find indices of events
  tradedata(removalindextime,:) = [];  %remove entries from data
  
  % generate time bars from tick data
  bars = Resampling003aa(barsize,tradedata);
  
%   %convert dates to numbers: Not needed!
%   timebars = zeros(size(bars));
%   timebars(:,1) = datenum(bars(:,1));
%   %dimensions = [datetime, VWAP, vol, Po, Ph, Pl, Pc]
%   timebars(:,2:7) = cell2mat(bars(:,2:7)); 
%   
  %[P_vwap,Vol,P_c,P_h,P_l,P_o,DateTime]
 
%   % Extract relevant variables: Open, Close, High, Low, Date-Time
%   DateTime = timebars(:,1);
%   P_vwap = cell2mat(timebars(:,2));
%   Vol = cell2mat(timebars(:,3));
%   P_o = cell2mat(timebars(:,4));
%   P_c = cell2mat(timebars(:,7));
%   P_h = cell2mat(timebars(:,5));
%   P_l = cell2mat(timebars(:,6));
  
  % write data to excel file
  xlswrite(char(strcat(writepath,datafilenames{l}(1:3),'-6M-ticks-timebars.xls')),bars);
  
  clearvars rawdatacell appendrawcell rawdata bars hr mn sc tradedat removalindextime
end

%% STEFI
%fields in raw data:
%index,times,type,value,size,condcode

[~,~,rawSTEFI] = xlsread('C:\Users\nicjm\Documents\MATLAB\Masters\Data\Bloomberg\tick\STEFI\STEFI-Index-6M.csv');
STEFIcell = regexp(rawSTEFI,',','split');
STEFI = vertcat(STEFIcell{:});
STEFI(:,[1,3,5]) = [];
STEFI(:,1) = cellstr(datestr(datenum(datestr(STEFI(:,1))),'yyyy-mm-dd')); %,'InputFormat','yyyy/mm/dd HH:MM');
  xlswrite(char(strcat(writepath,'\STEFI\STEFI-Index-6M.xlsx')),STEFI);
  
%% Naspers
rawdata = readtable(datafilenames{21}); %load csv
f = table2cell(rawdata); %convert to cell array from table
f1 = f(:,[2,4,5]); %extract needed cols only
clearvars f rawdata
for i = 1:size(f1,1)
    df(i) = datenum(f1{i,1});
end

npn = cell(size(f1));
%npn(:,1) = num2cell(df);
%npn(:,1) = df;
ds = cellstr(datestr(datetime(datevec(df(:)))+hours(2),'yyyy-mm-dd HH:MM:SS'));
npn(:,1) = ds;
npn(:,2:3) = f1(:,2:3);
clearvars f1 df

[~,~,~,hr,mn,sc] = datevec(ds); %find hour, minute and second
  %[tradedatetime] = datenum(yr,mt,dy,hr,mn,sc);  %create vector of dates and times 
    
 % appendrawcell(:,1) = cellstr(datestr(datetime(appendrawcell(:,1))+hours(2),'yyyy-mm-dd HH:MM:SS'));

% Remove all events between before 09h10 and after 16h50
removalindextime(1:size(((hr<9) | (hr==9&mn<=10) | (hr==16&mn>30&sc>=0) | (hr==17)),1)) = ((hr<9) | (hr==9&mn<=10) | (hr==16&mn>30&sc>=0) | (hr==17));
npn(removalindextime,:) = []; 
%npn(:,1) = datestr(str2double(npn(:,1))
bars = Resampling003aa(5,npn);

% write data to excel file
xlswrite(char(strcat(writepath,datafilenames{21}(1:3),'-6M-ticks-timebars.xls')),bars);
  
