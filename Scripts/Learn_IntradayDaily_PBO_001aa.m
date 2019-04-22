%% Masters: Online Learning for intraday data
%
% Author: N.J. Murphy
%
% Person Version:
% Masters-scripts-
%
% 1 Problem Specification: 
%
% 2 Data Specification: 6 months worth of transaction data for 20 stocks listed
%   on the JSE Top 40. The specified period is from 01 January 2018 to 30 June
%   2018. 
%
% 3 Configuration Control:
%        userpath/MATLAB/Masters
%        userpath/MATLAB/Masters       
%        userpath/MATLAB/Masters/Scripts    
%        userpath/MATLAB/Masters/Functions 
%        userpath/MATLAB/Masters/Data
%        userpath/MATLAB/Masters/html        
%
% 4 Version Control: No current version control
%
% 5 References: None
%
% 6 Current Situation: 
% 7 Future Situation: no future situation
%
% Uses: This script works with learning_IntradayDaily005. Does analysis for 
% both TC and no TC case.   

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Data Description %%

% 1  RIC               - Reuters instrument code.
%                              A Reuters instrument code, or RIC, is a ticker
%                              like code used by Thomson Reuters to identify
%                              financial instruments and indices.
% 
% 2   DateL            - Date of transaction.
% 3   TimeL            - Time of transaction.
% 4   DateTimeL        - Date and time of transaction.
% 5   Type             - Type of transaction: * Auction
%                                             * Quote
%                                             * Trade
% 6   Price            - Price of transaction for auction and trade
% 7   Volume           - Volume of tramsaction for auction and trade
% 8   MarketVWAP       - Market Volumes Weighted average prices
%                        = sum(shares bought * share price)/(Total shares bought)
% 9   L1BidPrice       - The price a buyer is willing to pay
% 10  L1BidSize        - The size of order a buyer wishes to achieve
% 11  L1AskPrice       - The price the seller is willing to offer
% 12  L1AskSize        - The size of order the seller wishes to achieve
%

%%%%%%%%%%%%%%%%%%%
%% 2. Data Cases %%

%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Clear Workspace %%
% This section prepares the workspace for the implementation of the script.
clc;            % clear command window
clear;
format long g;  % formating for output on comand window
format compact; % formating for output on comand window

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
addpath(fullfile(userpathstr,projectpath,'Functions','StatArb'));
addpath(fullfile(userpathstr,projectpath,'Functions','BCRP'));
addpath(fullfile(userpathstr,projectpath,'Functions','Indicators'));
addpath(fullfile(userpathstr,projectpath,'Scripts'));
addpath(fullfile(userpathstr,projectpath,'Data\Bloomberg\tick\Processed2'));

% Paths for strategies
stratpath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Functions\Strategies';
stratdirec = dir(fullfile(stratpath));
stratfilenames = {stratdirec(~[stratdirec.isdir]).name}';

% Paths for data
%Intraday:
datapath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Data\Bloomberg\tick\Processed2';
datadirec = dir(fullfile(datapath, '*.xls'));
datafilenames = {datadirec(~[datadirec.isdir]).name}';
%Daily:
datapathdaily = 'C:\Users\nicjm\Documents\MATLAB\Masters\Data\Bloomberg\Daily\6M';

%STEFI:
STEFIpath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Data\Bloomberg\tick\Processed\STEFI\STEFI-Index-6M.xlsx';
%datadirecdaily = dir(fullfile(datapathdaily));
%datafilenamesdaily = {datadirecdaily(~[datadirecdaily.isdir]).name}';

%%%%%%%%%%%%%%%%%%
%% 5. Load Data %%
% Load data for opening and closing prices, high and low prices, volume
% traded and dates.

nostrats = input('Number of strategies = ');
noofportfoliostocks = input('Number of stocks to be used in portfolio construction = ');
noofdays = input('Number of days = '); %the max number of days is 121 (123 if not worrying about liquidity calculation)
noofstocks = length(datafilenames);
dayperiodinterval = 4;%no. of days to use to pick liquid stocks
periodinterval = 88*dayperiodinterval; %no. of time bars to lookback for finding liquid stocks: 10 days worth
weighttype = 'volloading';

%extract data from .csv file and create Data matrix to be sent to learning
%class
[r,~] = size(xlsread(datafilenames{1}));%,range)); %get number of time bars for entire period
Data = nan(noofstocks,5,r);
[~,Datetimes_intra,~] = xlsread(datafilenames{1});%,strcat('A',num2str(periodinterval+1),':A);%,'A1:A1500');

for l = 1:length(datafilenames) 
  % Step 1: load processed data for entire period: dimensions = [datetime, VWAP, vol, Po, Ph, Pl, Pc]
  [~, ~, timebars] = xlsread(datafilenames{l});%,range); %load csv
  
  % deal with missing data
  [rowmiss,colmiss,~] = find(cell2mat(timebars(:,2:end))==0);
  timebars(rowmiss,colmiss+1) = {NaN};
  %fill in empty rows (time bars during which no transactions took place)
  barsnew = fillmissing(cell2mat(timebars(:,2:7)),'previous');
  barsnew = fillmissing(cell2mat(timebars(:,2:7)),'nearest'); %for missing values at time 1
  
  % Step 2: Create Data matrix- [P_c, P_o, P_h, P_l, P_v]
  Data(l,:,:) = [barsnew(:,6)'; barsnew(:,3)'; barsnew(:,4)'; barsnew(:,5)'; barsnew(:,2)'];
  clearvars barsnew timebars
end
  
%% Liquidity
%rank noofstocks most liquid stocks from the Top 40
Vol = reshape(Data(:,end,:),size(Data,1),size(Data,3))'; 

%get average daily volume for last 'periodinterval' periods
average5minvolperiodinterval = sum(Vol(1:periodinterval,:))/periodinterval;
[~,liquidstocks] = sort(average5minvolperiodinterval,'descend');

%define new data matrix of top 'noofportfoliostocks' most liquid stocks
DataIntraday = Data(:,:,periodinterval+1:end);
Datetimes_intra = Datetimes_intra(periodinterval+1:end,1);

% Get Bloomberg codes for stocks
stocks = extractfield(datadirec,'name');
stocknames = cell(1,noofstocks);
for i=1:noofstocks
   stocknames{i} = stocks{i}(1:3);
end

% Compute the sector clusters (3 sectors + trivial)
sectors = sectors_intraday(stocknames(liquidstocks(1:noofportfoliostocks))); %organise the 3 sector clusters
sectors{4} = ones(1,noofportfoliostocks); % trivial cluster (all stocks)

%check to see if there are any empty clusters
for i = 1:4
    if any(sectors{i})==0
      sectors{i} = nan;
    end
end
sectors(cellfun(@(sectors) any(isnan(sectors)),sectors)) = []; %remove nan's

clearvars average5minvolperiodinterval stocks Vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Daily Data
%load daily OHLCV data from file
[~,~,P_c] = xlsread(strcat(datapathdaily,'\JSEClosing'));%,'B2:AD125');
[~,~,P_o] = xlsread(strcat(datapathdaily,'\JSEOpening'));%,'B2:AD125');
[~,~,P_v] = xlsread(strcat(datapathdaily,'\JSEVol'));%,'B2:AD125');
[~,~,P_h] = xlsread(strcat(datapathdaily,'\JSEHigh'));%,'B2:AD125');
[~,~,P_l] = xlsread(strcat(datapathdaily,'\JSELow'));%,'B2:AD125');

%Remove stock names and dates from cell array
datesdaily = P_c(2+dayperiodinterval:end,1);
stocknamesdaily = P_c(1,2+dayperiodinterval:end);
P_c(1,:) = [];
P_o(1,:) = [];
P_h(1,:) = [];
P_l(1,:) = [];
P_v(1,:) = [];

P_c(:,1) = [];
P_o(:,1) = [];
P_h(:,1) = [];
P_l(:,1) = [];
P_v(:,1) = [];

% Aggreagte data into matrix
Datadaily = nan(noofstocks,5,noofdays);
for i = 1:noofstocks
  Datadaily(i,:,:) = [cell2mat(P_c(dayperiodinterval+1:noofdays+dayperiodinterval,i)'); cell2mat(P_o(dayperiodinterval+1:noofdays+dayperiodinterval,i)');...
      cell2mat(P_h(dayperiodinterval+1:noofdays+dayperiodinterval,i)'); cell2mat(P_l(dayperiodinterval+1:noofdays+dayperiodinterval,i)');...
      cell2mat(P_v(dayperiodinterval+1:noofdays+dayperiodinterval,i)')];
end

% remove the data up to noofdays from Intraday set
[~,ia1,intradaycutoff] = intersect(string(datesdaily{noofdays+1}),datestr(datevec(Datetimes_intra),'yyyy-mm-dd'));
DataIntraday = DataIntraday(:,:,1:intradaycutoff-1);

% find indices for unique days, start of days and end of days
[yr,mt,dy,~,~,~] = datevec(Datetimes_intra(1:intradaycutoff));
tradedate = datenum(yr,mt,dy);
[~,uniquedaytr] = unique(tradedate,'first');  %find the indices for each day 

Datetimes_intra = Datetimes_intra(1:intradaycutoff-1);


clearvars l i hr dy mn mt P_o P_v P_h P_l rowmiss colmiss tradedate yr r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load STEFI
[~,~,STEFIraw] = xlsread(STEFIpath);
[~,ia,~] = intersect(STEFIraw(:,1), datesdaily(1:noofdays)); %pick out the 
STEFI = str2double(STEFIraw(ia,2));

%% Prepare strategies
% Specify k and ell values to be used as parameters and the number of
% strategies
%usual one:
ell = 3:6:54;
k = 6:6:64;

% Specify the number of strategies  
strats = cell(nostrats,3);
strats(1:nostrats,1) = stratfilenames(1:nostrats,1);

% Determine the names of input parameters to see how many
% free parameters are re.quired
for i = 1:nostrats
    strategy = strats{i}(1:end-2);
    strategyfn = str2func(strategy);
    fid = fopen(strcat(strategy,'.m')); %read the 1st line of the cuurent strategy function
    tline = fgetl(fid);              
    fclose(fid);
    Inputs = tline(strfind(tline, '(')+1:end-1); %this can just go into next line
    InputsVec = strsplit(Inputs,', '); %names of function input args 
    strats{i,2} = any(cell2mat(strfind(InputsVec ,'n2'))); %parameter condition to check if n2 is an input
    strats{i,3} = InputsVec;
end

% Sort the strategies in an order such that the strategies with 2 params
% come first
[~,order] = sort(cell2mat(strats(:,2)),'descend');
strats = strats(order,:);

clearvars stratfilenames strategyfn strategy tline order InputsVec Inputs ...
    fid i ia range STEFIraw stratdirec stratpath STEFIpath
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%-------------------------%%
%% 7. Monte-Carlo Learning %%
%%-------------------------%%
inds = [1:3:45,45];
%inds = [1:3:10];


for i = 1:length(inds)-1
    % compute indices for data piece to cut out
    ind_leave_out_daily = inds(i):(inds(i+1)-1);
    ind_leave_out_intra = 88*(inds(i)-1)+1:88*(inds(i+1)-1);
    ind_jump = ind_leave_out_daily(end)-2; %daily jump is needed only
    
    %remove data from both daily and intraday sets
    DataIntraday_temp = DataIntraday;
    DataIntraday_temp(:,:,ind_leave_out_intra) = []; 
    DataDaily_temp = Datadaily;
    DataDaily_temp(:,:,ind_leave_out_daily) = [];
    
    % Prepare data for learning by calling 'learning' function from 'learning' class
    p = learning_IntradayDaily_PBO(DataIntraday_temp,DataDaily_temp,k,ell,strats,STEFI,periodinterval,stocknames,sectors,noofportfoliostocks,liquidstocks,weighttype,uniquedaytr,ind_jump);

    % Implement the learning procedure by calling 'offline' function from the 'learning' class
    p1 = offline_Intraday(p);
    
    M_pnl(:,i) = p1.PnL;
    M_S(:,i) = p1.S;
end


%create folder to save plota and workspace
%savepath = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC'));
%savefolder = char(strcat('New',num2str(noofportfoliostocks),'stocks','_',num2str(nostrats),'strats','_days',num2str(noofdays)));
%mkdir(strcat(savepath,'\',savefolder))
xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\M_Tim_pnl','.xlsx'),M_pnl) 
xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\M_Tim_S','.xlsx'),M_S) 



% %%%%%%%%%%%%%%%%%%%%%
% %% 8. Plot Results %%
% % get the names of the strategies used
% strat = cell(1,nostrats); 
% for j = 1:size(strats,1)
%  strat{1,j} = strats{j,1}(1:end-2);
% end
% 
% clearvars i j ia1 p 
% % first date of active trading
% tmin0 = uniquedaytr(2);
% 
% % Convert dates to date time vector
% Datetimevec = datetime(datesdaily(1:noofdays));
% x = linspace(Datetimevec(1),Datetimevec(end),length(p1.S(tmin0:end)));
% 
% % Specify plot information
% plotVar = {'p1.S(tmin0:end,:)','cumsum(p1.PnL(tmin0:end,:))','p1.SH_fused(tmin0:end,:)','p1.b(tmin0:end,:)','p1.SHave_fused(tmin0:end,:)'...
%     ,'cumsum(p1.PnL_experts(tmin0:end,:))','p1.S_TC(tmin0:end,:)','cumsum(p1.PnL_TC(tmin0:end,:))'};
% plotTitle = {'Overall Portfolio Wealth ({\bf{S}}) and Profits and Losses ({\bf{PL}})','Cumulative Profits and Losses','Relative Population Wealth of Experts','Overall Portfolio Controls','Relative Population Wealth of Strategies',...
%     'Cumulative Profits and Losses of experts','Overall Portfolio Wealth ({\bf{S}}) and Profits and Losses ({\bf{PL}})','Cumulative Profits and Losses'};
% plotX = {'Time'};
% plotY = {'Wealth ({\bf{S}})','{\bf{PL}}','Expert Wealth ({\bf{Sh}_{t}})','Weight ({\bf{b}})','Mean Expert Wealth per Strategy','{\bf{PL}}','Wealth ({\bf{S}})','{\bf{PL}}'};
% Savename = {'S','PnL','SH','b','SHave','PnL_experts','S_TC','PnL_TC'};
% 
% %get proper strategy names
% strategy_names = ["EMA X over";"Ichimoku Kijun Sen";"MACD";"Moving Ave X over";"ACC";"Bollinger";"Fast Stochastic";"MARSI";"MOM";"Anti-BCRP";"Anticor";"Online BCRP";"PROC";"RSI";"SAR";"Slow Stochastic";"Williams %R"];
% 
% for i = 1:length(plotVar)-1
%     
%     if i==1 || i==7
%         j = i;
%         hf1=figure;
%         var = eval(plotVar{j});
%         plot(x,var);
%         datetick('x','dd-mmm-yyyy','keepticks','keeplimits')
% 
%         hold on 
%         set(gca,'XTickLabelRotation',50)  
%         title(plotTitle{j});
%         ylabel(plotY{j});
%         %xlabel(plotX{1});
%         if j == 7
%           ylim([0,7])
%         end
%         axes('parent',hf1,'position',[0.2 0.6 0.3 0.3]);
%         j=j+1;
%         var = eval(plotVar{j});
%         plot(x,var);
%         %grid on;
%         %set(gca,'xticklabel',{[]})
%         datetick('x','dd-mmm-yyyy','keepticks','keeplimits')
%         
%         hold on 
%         set(gca,'XTickLabelRotation',50,'FontSize',6)  
%         
%         %title(plotTitle{i});
%         ylabel(plotY{j});
%         if j == 8
%           ylim([-0.05,2])
%         else
%            ylim([0,2])
%         end
%         %xlabel(plotX{1});
%         hold off
%             
%           %%%%%%%%%%%%%%%
%           % save figure %
%           %%%%%%%%%%%%%%%
%            savetoname = char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_',Savename{i}));
%            fig = gcf;
%            fig.PaperUnits = 'centimeters';
%            %fig.PaperPosition = [0 0 19 12]; %proposal size
%            fig.PaperPosition = [0 0 12 12];
%            print(savetoname,'-dpng','-r300')
%     else
%         figure;
%         var = eval(plotVar{i});
%         plot(x,var);
%         datetick('x','dd-mmm-yyyy','keepticks','keeplimits')
% 
%         hold on 
%         set(gca,'XTickLabelRotation',50)  
%         title(plotTitle{i});
%         ylabel(plotY{i});
%         %xlabel(plotX{1});
%         hold off
%             
%           %%%%%%%%%%%%%%%
%           % save figure %
%           %%%%%%%%%%%%%%%
%            savetoname = char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_',Savename{i}));
%            fig = gcf;
%            fig.PaperUnits = 'centimeters';
%            %fig.PaperPosition = [0 0 19 12]; %proposal size
%            fig.PaperPosition = [0 0 12 12];
%            print(savetoname,'-dpng','-r300')
%     end
%     
%     if i == 8 || i == 7 
%             figure;
%         var = eval(plotVar{i});
%         plot(x,var);
%         datetick('x','dd-mmm-yyyy','keepticks','keeplimits')
% 
%         hold on 
%         set(gca,'XTickLabelRotation',50)  
%         title(plotTitle{i});
%         ylabel(plotY{i});
%         %xlabel(plotX{1});
%         hold off
%             
%           %%%%%%%%%%%%%%%
%           % save figure %
%           %%%%%%%%%%%%%%%
%            savetoname = char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_',Savename{i},'_disser'));
%            fig = gcf;
%            fig.PaperUnits = 'centimeters';
%            %fig.PaperPosition = [0 0 19 12]; %proposal size
%            fig.PaperPosition = [0 0 12 12];
%            print(savetoname,'-dpng','-r300')
%      end
% 
%  
% end
% clearvars fig
%  
% % prepare variables for saving
% S = p1.S(tmin0:end,:);
% SH_fused = p1.SH_fused(tmin0:end,:);
% SHave_fused = p1.SHave_fused(tmin0:end,:);
% b = p1.b(tmin0:end,:);
% PnL_experts =  p1.PnL_experts(tmin0:end,:);
% PnL = p1.PnL(tmin0:end,:);
% PnL_TC = p1.PnL_TC(tmin0:end,:);
% TC_vec = p1.TC_vec(tmin0:end,:);
% 
% %call expert parameter function to get parameters of each expert
% parameters = ExpertParameters(4,ell,k,strats);
% 
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_stockNames.xlsx')),stocknames(liquidstocks(1:noofportfoliostocks)));
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_PnLexperts.xlsx')),p1.PnL_experts);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_PnL.xlsx')),p1.PnL);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_PnL_TC.xlsx')),p1.PnL_TC);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_SH_fused.xlsx')),p1.SH_fused);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_S.xlsx')),p1.S);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_S_TC.xlsx')),p1.S_TC);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_SHave_fused.xlsx')),p1.SHave_fused);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_STEFI.xlsx')),STEFI);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_parameters.xlsx')),parameters);
% xlswrite(char(strcat(userpathstr,'\Masters\Plots\IntradayDaily2one\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_TC.xlsx')),TC_vec);

%-------------------------------------------------------------------------------------