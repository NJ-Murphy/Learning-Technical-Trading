%% Masters: Online Learning for daily data
%
% Author: N.J. Murphy
%
% 1 Problem Specification: 
%
% 2 Data Specification: Set of OHLCV data values for 42 stocks on the JSE Top 40
%   over the period 01-01-2005 to 29-04-2016
%
% 3 Configuration Control:
%        userpath/MATLAB/Masters
%        userpath/MATLAB/Masters/Scripts    
%        userpath/MATLAB/Masters/Functions 
%        userpath/MATLAB/Masters/Data    
%
% 4 Version Control: No current version control
%
% 5 References: None
%
% 6 Current Situation: 
% 7 Future Situation: no future situation
%
% Uses: follows from 013aa and works with 019. changes 013 noofdays
% variable and rather looks back on data before noofdays to get liquid
% stocks.
%
%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Clear Workspace %%
% This section prepares the workspace for the implementation of the script.
clc;            % clear command window
format long g;  % formating for output on comand window
format compact; % formating for output on comand window
clear % clear workspace

%%%%%%%%%%%%%%%%%%%
%% 2. Path Setup %%
% Set the project paths
[projectPath] = pwd;

% Add the project path so that the script sees the necessary functions
addpath(projectPath);
userpathstr = userpath;
userpathstr = userpathstr(~ismember(userpathstr,';'));
projectpath = 'Masters';   %path for personal computer
addpath(fullfile(userpathstr,projectpath,'Functions'));
addpath(fullfile(userpathstr,projectpath,'Functions\StatArb'));
addpath(fullfile(userpathstr,projectpath,'Functions\BCRP'));
addpath(fullfile(userpathstr,projectpath,'Functions','Strategies'));
addpath(fullfile(userpathstr,projectpath,'Functions','Indicators'));
addpath(fullfile(userpathstr,projectpath,'Scripts'));
addpath(fullfile(userpathstr,projectpath,'Data\TRTH\EOD'));

% # Filepath to input data
[datafilepath] = 'C:\Users\nicjm\Documents\MATLAB\Masters\Data\TRTH\EOD';
datadirec = dir(fullfile(datafilepath));
datafilenames = {datadirec(~[datadirec.isdir]).name}';

% Paths for strategies
stratpath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Functions\Strategies';
stratdirec = dir(fullfile(stratpath));
stratfilenames = {stratdirec(~[stratdirec.isdir]).name}';

%%%%%%%%%%%%%%%%%%
%% 3. Load Data %%
% Load data for opening and closing prices, high and low prices, volume
% traded and dates.

% Load stocks: 
noofdays = input('Number of days = '); %max is 4125
nostrats = input('Number of strategies = ');
noofportfoliostocks = input('Number of stocks to be used in portfolio construction = ');
dayinterval = 375;
weighttype = 'volloading'; %'invvolloading'; 

% Extract the OHLCV stock data for all stocks for the required period
% 'noofdays'. Here we go from most recent (4139) back a number of
[L,~] = size(xlsread('JSEClosing',strcat('C4139:AS',num2str(4138-noofdays))));
while L ~= noofdays
  noofdays = noofdays + 1;
  [L,~] = size(xlsread('JSEClosing',strcat('C4139:AS',num2str(4138-noofdays))));
end  
P_c(1:noofdays,:) = xlsread('JSEClosing',strcat('C4139:AS',num2str(4138-noofdays)));
P_o(1:noofdays,:) = xlsread('JSEOpening',strcat('C4139:AS',num2str(4138-noofdays)));
P_v(1:noofdays,:) = xlsread('JSEVol',strcat('C4139:AS',num2str(4138-noofdays)));
P_l(1:noofdays,:) = xlsread('JSELow',strcat('C4139:AS',num2str(4138-noofdays)));
P_h(1:noofdays,:) = xlsread('JSEHigh',strcat('C4139:AS',num2str(4138-noofdays)));

% error('Please choose different number of days as last day falls on weekend/public holiday')

[~,stocknames,] = xlsread('JSEClosing','C2:AS2');

% Load dates
[~,~,Dates] = xlsread('JSEClosing',strcat('A4139:A',num2str(4138-noofdays+2)));

% Load STEFI-> 2864 is the last date in the OHLCV daily data
[STEFI,~,~] = xlsread('PT-TAA-JSE-Daily-1994-2017',2,strcat('C253:C',num2str(253+noofdays-1)));
STEFI = fliplr(STEFI')'; 

clearvars stratpath
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Process and Prepare the Data for learning %%
% Here, we prepare the data for input into the learning algorithm

% Remove any stocks for which there are more than 80% of the data point
% missing
[r,c] = size(P_v); 

for i = 1:c
    if sum(isnan(P_c(:,i)))/r > 0.50  %If more than 50% of the column is NAN then remove this column
      y(i,1) = i;
    else
      y(i,1) = 0;  
    end
end
z = y(y~=0);
P_c(:,z(1:end)) = [];
P_o(:,z(1:end)) = [];
P_h(:,z(1:end)) = [];
P_l(:,z(1:end)) = [];
P_v(:,z(1:end)) = [];

% Compute number of stocks
[~,noofstocks] = size(P_c);
stocknames(z(1:end)) = [];

% Get RIC for stocks
RIC = cell(1,noofstocks);
for i=1:noofstocks
   RIC{i} = stocknames{1,i}(1:end-3);
end

% Aggregate opening and closing prices, high and low prices, volume
% traded and dates into one single array (Data = [Stocks,Prices,Time]).
Data = nan(noofstocks,5,r);
for i = 1:noofstocks
  Data(i,:,:) = [P_c(:,i)'; P_o(:,i)'; P_h(:,i)'; P_l(:,i)'; P_v(:,i)'];
end

clearvars P_o P_h P_l P_v y i c r stratdirec projectpath projectPath 
% Remove NaNs from data (weekends and public holidays)
Dates = Dates(1:size(Data,3)); %remove dates which 'Data' did not import due to empty cells
Dates(any(any(isnan(Data),1),2)) = []; %remove dates of stocks which didnt trade the whole period, also weekends
STEFI(reshape(any(any(isnan(Data),1),2),1,length(P_c))) = [];
Data(:,:,any(any(isnan(Data),1),2)) = [];

% check to see STEFI and Data have the same number of entries:
if size(STEFI,1)~=size(Data,3)
   error('Not enough data for STEFI, please change the number of days')
end

%replace NaNs in STEFI with interpolated values - missing data?
for i = 1:length(STEFI)
   if isnan(STEFI(i)) == 1
     STEFI(i) = ((i-(i-1))/(i+1-(i-1)))*STEFI(i+1)+ ((i+1-i)/(i+1-(i-1)))*STEFI(i-1);
   end
end

% Specify k and ell values to be used as look-back parameters
ell = 4:6:54;
k = 6:6:64;

% capture the names of strategies  
strats = cell(nostrats,3);
strats(1:nostrats,1) = stratfilenames(1:nostrats,1);

% Determine the names of input parameters to see how many
% free parameters are recquired
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

%%%%%%%%%%%%%%%%%%
%% 5. Liquidity %%
%find 'noofportfoliostocks' most liquid stocks from the Top 40 to be traded
%for the first 120 days

% get volume from dayinterval days before active trading starts - check
% that theres sufficent data
if 4138-noofdays-dayinterval-1 < 0 
    error('choose fewer ''noofdays'' argument')
end
Volume_liquid = xlsread('JSEVol',strcat('C',num2str(4138-noofdays-1-dayinterval),':AS',num2str(4138-noofdays)));
Volume_liquid(:,z(1:end)) = [];

% ADV
averagedailyvol = mean(Volume_liquid,1,'omitnan');
averagedailyvol(isnan(averagedailyvol)) = 0;
[~,liquidstocks] = sort(averagedailyvol,'descend');

% do the sector clusters (3 sectors + trivial)
sectors = sectors_daily(stocknames(liquidstocks(1:noofportfoliostocks))); %organise the 3 sector clusters
sectors{4} = ones(1,noofportfoliostocks); % trivial cluster (all stocks)

%check to see if there are any empty clusters
for i = 1:4
    if any(sectors{i})==0
      sectors{i} = nan;
    end
end
sectors(cellfun(@(sectors) any(isnan(sectors)),sectors)) = []; %remove nan's

clearvars stratfilenames z strategyfn strategy tline order InputsVec Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first date of active trading: Below the first dayinterval days are used
% for finding liquid stocks. The last tmin0 of those days are then gathered
% as lookback data which are sufficient for the longest lookback needed by
% any of the strategies. 

% Prepare data for learning by calling 'learning' function from 'learning' class
p = learning_daily(Data,k,ell,strats,STEFI,dayinterval,stocknames,sectors,noofportfoliostocks,liquidstocks,weighttype);

%%%%%%%%%%%%%%%%%
%% 6. Learning %%
% Input the processed data into the learning algorithm 

% Implement the learning procedure by calling 'offline' function from the 'learning' class
p1 = offline(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Best stock and uniform BAH
% Best Stock
% Get closing prices and log returns/price relatives
% P_c = reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3)); %(stocks*time)
% xrel = P_c(:,2:end)./P_c(:,1:end-1);
% 
% % find best stock 
% [~,beststockind] = max(prod(xrel,2));
% 
% % compute wealth
% S_beststock(1) = 1;
% S_beststock = cumprod([S_beststock(1) ,xrel(beststockind,:)]);
% 
% % BAH
% b_BH = repmat(1/noofstocks,1,noofstocks);
% bb = repmat(b_BH,size(xrel,2),1);
% %S_uBAH(1) = 1;
% S_uBAH = [1,cumprod(sum(bb'.*xrel),1)];

%%%%%%%%%%%%%
%% 8. BCRP %%
% simulation of best offline CRP

%create folder to save plota and workspace
savepath = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC'));
savefolder = char(strcat('New',num2str(noofportfoliostocks),'stocks','_',num2str(nostrats),'strats','_days',num2str(noofdays)));
mkdir(strcat(savepath,'\',savefolder))

Pc = reshape(Data(liquidstocks(1:noofportfoliostocks),1,1:end),size(Data(liquidstocks(1:noofportfoliostocks),1,:),1),size(Data(liquidstocks(1:noofportfoliostocks),1,1:end),3)); %closing prices time t-1 to t          
ret_all = [transpose(Pc(:,2:end)./Pc(:,1:end-1)),STEFI(2:end,1)./STEFI(1:end-1,1)];

%BCRP_opt(noofportfoliostocks+1,ret_all)
%R test: xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\returns.xlsx'),ret_all) 
[BCRP_perf,best_port] = BCRP_daily(ret_all,5000);
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',num2str(noofdays),'BCRP','.xlsx'),BCRP_perf,1) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',num2str(noofdays),'BCRP','.xlsx'),best_port,2) 


%%%%%%%%%%%%%%%%%%%%%
%% 9. Plot Results %%

% get the names of the strategies used
strat = cell(1,nostrats); 
for j = 1:size(strats,1)
 strat{1,j} = strats{j,1}(1:end-2);
end

% Convert dates to date time vector
tmin0 = 1;
Datetimevec = datetime(Dates(tmin0:end,:));

% Specify plot information
plotVar = {'p1.S(tmin0:end,:)','cumsum(p1.PnL(tmin0:end,:))','p1.SH(tmin0:end,:)','p1.b(tmin0:end,:)','p1.SHave(tmin0:end,:)'...
    ,'cumsum(p1.PnL_experts(tmin0:end,:))','p1.S_TC(tmin0:end,:)','cumsum(p1.PnL_TC(tmin0:end,:))'};
plotTitle = {'Overall Portfolio Wealth ({\bf{S}}) and Profits and Losses ({\bf{PL}})','Cumulative Profits and Losses','Relative Population Wealth of Experts','Overall Portfolio Controls','Relative Population Wealth of Strategies',...
    'Cumulative Profits and Losses of experts','Overall Portfolio Wealth ({\bf{S}}) and Profits and Losses ({\bf{PL}})','Cumulative Profits and Losses'};
plotX = {'Time'};
plotY = {'Wealth ({\bf{S}})','{\bf{PL}}','Expert Wealth ({\bf{Sh}_{t}})','Weight ({\bf{b}})','Mean Expert Wealth per Strategy','{\bf{PL}}','Wealth ({\bf{S}})','{\bf{PL}}'};
Savename = {'S','PnL','SH','b','SHave','PnL_experts','S_TC','PnL_TC'};

%get proper strategy names
%strategy_names = cell(size(strats,1));
strategy_names = ["EMA X over";"Ichimoku Kijun Sen";"MACD";"Moving Ave X over";"ACC";"Bollinger";"Fast Stochastic";"MARSI";"MOM";"Anti-BCRP";"Anticor";"Online BCRP";"PROC";"RSI";"SAR";"Slow Stochastic";"Williams %R"];

for i = 1:length(plotVar)-1
    
   if i==1 || i==7
        j = i;
        hf1=figure;
        var = eval(plotVar{j});
        plot(Datetimevec,var(2:end,:));
        datetick('x','dd-mmm-yyyy','keepticks','keeplimits')

        hold on 
        set(gca,'XTickLabelRotation',50)  
        title(plotTitle{j});
        ylabel(plotY{j});

        if i == 1
          plot(Datetimevec,BCRP_perf)
          legend('OLA',' BCRP','Location','southeast')
        end

        %xlabel(plotX{1});
        if j == 7 || j==1
          ylim([0,6])
        end
        axes('parent',hf1,'position',[0.2 0.6 0.3 0.3]);
        j=j+1;
        var = eval(plotVar{j});
        plot(Datetimevec,var(2:end,:));
        %grid on;
        %set(gca,'xticklabel',{[]})
        datetick('x','dd-mmm-yyyy','keepticks','keeplimits')

        hold on 
        set(gca,'XTickLabelRotation',50,'FontSize',6)  

        %title(plotTitle{i});
        ylabel(plotY{j});
        if j == 8
          ylim([-0.05,2])
        else
           ylim([0,2])
        end
        %xlabel(plotX{1});
        hold off

        %%%%%%%%%%%%%%%
        % save figure %
        %%%%%%%%%%%%%%%
        savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_',Savename{i}));
        fig = gcf;
        fig.PaperUnits = 'centimeters';
        %fig.PaperPosition = [0 0 19 12]; %proposal size
        fig.PaperPosition = [0 0 12 12];
        print(savetoname,'-dpng','-r300')
   else
        figure;
        var = eval(plotVar{i});
        plot(Datetimevec,var(2:end,:));
        datetick('x','dd-mmm-yyyy','keepticks','keeplimits')

        hold on 
        set(gca,'XTickLabelRotation',50)  
        title(plotTitle{i});
        ylabel(plotY{i});
        %xlabel(plotX{1});
        ylim([-0.1,1.4])
        hold off
        
        %%%%%%%%%%%%%%%
        % save figure %
        %%%%%%%%%%%%%%%
        savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_',Savename{i}));
        fig = gcf;
        fig.PaperUnits = 'centimeters';
        %fig.PaperPosition = [0 0 19 12]; %proposal size
        fig.PaperPosition = [0 0 12 12];
        print(savetoname,'-dpng','-r300')
   end
end

%stock names used in learning
stocks = stocknames(liquidstocks(1:noofportfoliostocks));
clearvars var strat savetoname Savename savepath RIC plotVar plotX plotY plotTitle...
          fid bb i j p

% prepare variables for saving
S = p1.S(tmin0:end,:);
S_TC = p1.S_TC(tmin0:end,:);
SH = p1.SH(tmin0:end,:);
SHave = p1.SHave(tmin0:end,:);
b = p1.b(tmin0:end,:);
PnL_experts =  p1.PnL_experts(tmin0:end,:);
PnL = p1.PnL(tmin0:end,:);
PnL_TC = p1.PnL_TC(tmin0:end,:);

%call expert parameter function to get parameters of each expert
parameters = ExpertParameters(4,ell,k,strats);

% save results:
% R folder
% mkdir(strcat('C:\Users\nicjm\Documents\R\Masters\',savefolder))
% xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\',savefolder,'\SH','.xlsx'),SH) 
% xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\',savefolder,'\S','.xlsx'),S) 
% xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\',savefolder,'\b','.xlsx'),b) 
% xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\',savefolder,'\params','.xlsx'),parameters) 
% xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\',savefolder,'\SHave','.xlsx'),SHave) 
% xlswrite(strcat('C:\Users\nicjm\Documents\R\Masters\',savefolder,'\Stefi','.xlsx'),STEFI(tmin0:end,:)) 

% Matlab folder
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'SH_TC','.xlsx'),SH) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'S','.xlsx'),S) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'S_TC','.xlsx'),S_TC) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'SHave_TC','.xlsx'),SHave) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'PnL_TC','.xlsx'),PnL_TC) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'PnL','.xlsx'),PnL) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'PnL_experts_TC','.xlsx'),PnL_experts) 
xlswrite(char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_stockNames.xlsx')),stocknames(liquidstocks(1:noofportfoliostocks)));
xlswrite(char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\stocks',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_TransCosts.xlsx')),p1.TC_vec);

%save workspace
%save(char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_days',string(noofdays))))


%% 10. Process SH and expertparams to be saved to file for use in Python

% Write to file
%xlswrite(strcat('C:\Users\nicjm\Documents\Python\VAEs\ExpertWealth','.xlsx'),SH) 
%xlswrite('C:\Users\nicjm\Documents\Python\VAEs\Params.xlsx',parameters) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11. Statistical Arbitrage Test no TC %%

%Datetimevec={'22-Oct-2009';'29-Apr-2016'}
%PnL = xlsread("C:/Users/nicjm/Documents/MATLAB/Masters/Plots/DailyDatawithTC/New15stocks_17strats_days3001/15_strats17_22-Oct-2009-29-Apr-2016_3001PnL.xlsx"); 
Statarb_S = cumsum(PnL(30:430));

% compute all required variable using stat arb function (CM)
[Min_t,t_c_95,p_val_95,p_loss,Min_tsim] = StatArb_cm(Statarb_S);

% plot and save histogram and probability of loss
figure;
bar(p_loss(1:20))
title('Probability of loss per period')
xlabel('n')
ylabel(strcat('{\bf Pr}','[loss after n periods]'),'Interpreter','tex')

savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_p_loss'));
fig = gcf;
fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 19 12]; %proposal size
fig.PaperPosition = [0 0 12 12];
print(savetoname,'-dpng','-r300')

% plot histogram of simulated Min-t's
figure;
histmint = histogram(Min_tsim,60);
xlabel('Min-$t$','Interpreter','latex')
ylabel('count')
title('{\bf Histogram of simulated Min-$t$''s from Monte Carlo}','Interpreter','latex')
line([t_c_95 t_c_95], ylim, 'Color', 'r')
%line([t_c_90 t_c_90], ylim, 'Color', 'r')
line([Min_t Min_t], ylim, 'Color', 'g')
hold on;
text(t_c_95+0.1, 260, '$t_{c}$','Interpreter','latex');
text(Min_t-0.6, 260, 'Min-$t$','Interpreter','latex');
hold off;

savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_Hist'));
fig = gcf;
fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 19 12]; %proposal size
fig.PaperPosition = [0 0 12 12];
print(savetoname,'-dpng','-r300')


% plot figures for pre-print (p-loss as inset)
h1=figure;
histogram(Min_tsim,60);
xlabel('Min-$t$','Interpreter','latex')
ylabel('count')
title('{\bf Histogram of simulated Min-$t$''s from Monte Carlo}','Interpreter','latex')
line([t_c_95 t_c_95], ylim, 'Color', 'r')
%line([t_c_90 t_c_90], ylim, 'Color', 'r')
line([Min_t Min_t], ylim, 'Color', 'g')
hold on
text(t_c_95+0.1, 260, '$t_{c}$','Interpreter','latex');
text(Min_t-0.6, 260, 'Min-$t$','Interpreter','latex');

axes('parent',h1,'position',[0.22 0.7 0.18 0.18]);
bar(p_loss(1:20))
%title('Probability of loss per period')
xlabel('n')
xlim([0 20])
set(gca,'FontSize',8)
ylabel(strcat('{\bf Pr}','[loss after n periods]'),'Interpreter','tex')
hold off

savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_Hist_pre_print'));
fig = gcf;
fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 19 12]; %proposal size
fig.PaperPosition = [0 0 12 12];
print(savetoname,'-dpng','-r300')


% save results in cell array and then to excel
stat_arb_out = cell(5,1);
stat_arb_out{1} = Min_t;
stat_arb_out{2} = t_c_95;
stat_arb_out{3} = p_val_95;
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'StatArb_results','.xlsx'),stat_arb_out,1) 

% save prob loss vector and min_t sim vector to excel doc
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'StatArb_results','.xlsx'),p_loss,2) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'StatArb_results','.xlsx'),Min_tsim,3) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 12. Statistical Arbitrage Test TC %%

%PnL_TC = xlsread("C:/Users/nicjm/Documents/MATLAB/Masters/Plots/DailyDatawithTC/New15stocks_17strats_days3001/15_strats17_22-Oct-2009-29-Apr-2016_3001PnL_TC.xlsx"); 
%S = xlsread("C:/Users/nicjm/Documents/MATLAB/Masters/StatArb/10_strats16_10-Dec-2007-29-Apr-2016_4100S.xlsx"); 
Statarb_S = cumsum(PnL_TC(30:430));

% compute all required variable using stat arb function (CM)
[Min_t,t_c_95,p_val_95,p_loss,Min_tsim] = StatArb_cm(Statarb_S);

% plot and save histogram and probability of loss
figure;
bar(p_loss)
title('Probability of loss per period')
xlabel('n')
ylabel(strcat('{\bf Pr}','[loss after n periods]'),'Interpreter','tex')

savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_p_loss_TC'));
fig = gcf;
fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 19 12]; %proposal size
fig.PaperPosition = [0 0 12 12];
print(savetoname,'-dpng','-r300')

% plot histogram of simulated Min-t's
figure;
histmint = histogram(Min_tsim,60);
xlabel('Min-$t$','Interpreter','latex')
ylabel('count')
title('{\bf Histogram of simulated Min-$t$''s from Monte Carlo}','Interpreter','latex')
line([t_c_95 t_c_95], ylim, 'Color', 'r')
%line([t_c_90 t_c_90], ylim, 'Color', 'r')
line([Min_t Min_t], ylim, 'Color', 'g')
hold on;
text(t_c_95+0.1, 270, '$t_{c}$','Interpreter','latex');
text(Min_t-0.7, 270, 'Min-$t$','Interpreter','latex');
hold off;

savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_Hist_TC'));
fig = gcf;
fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 19 12]; %proposal size
fig.PaperPosition = [0 0 12 12];
print(savetoname,'-dpng','-r300')


% plot figures for pre-print (p-loss as inset)
h1=figure;
histogram(Min_tsim,60);
xlabel('Min-$t$','Interpreter','latex')
ylabel('count')
title('{\bf Histogram of simulated Min-$t$''s from Monte Carlo}','Interpreter','latex')
line([t_c_95 t_c_95], ylim, 'Color', 'r')
%line([t_c_90 t_c_90], ylim, 'Color', 'r')
line([Min_t Min_t], ylim, 'Color', 'g')
hold on
text(t_c_95+0.1, 270, '$t_{c}$','Interpreter','latex');
text(Min_t-0.6, 270, 'Min-$t$','Interpreter','latex');

axes('parent',h1,'position',[0.22 0.7 0.18 0.18]);
bar(p_loss(1:400))
%title('Probability of loss per period')
xlabel('n')
xlim([0 400])
set(gca,'FontSize',8)
ylabel(strcat('{\bf Pr}','[loss after n periods]'),'Interpreter','tex')
hold off

savetoname = char(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'_Hist_pre_print_TC'));
fig = gcf;
fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 19 12]; %proposal size
fig.PaperPosition = [0 0 12 12];
print(savetoname,'-dpng','-r300')


% save results in cell array and then to excel
stat_arb_out = cell(5,1);
stat_arb_out{1} = Min_t;
stat_arb_out{2} = t_c_95;
stat_arb_out{3} = p_val_95;
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'StatArb_results_TC','.xlsx'),stat_arb_out,1) 

% save prob loss vector and min_t sim vector to excel doc
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'StatArb_results_TC','.xlsx'),p_loss,2) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'StatArb_results_TC','.xlsx'),Min_tsim,3) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 13. summary stats of expert and strtagey performance %%
[stats,stats_pnl]= expert_perf_stats(SH,PnL_experts,parameters,strategy_names);

xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'expert_stats_TC','.xlsx'),stats_pnl,1) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'expert_stats_TC','.xlsx'),stats,2) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 14. BCRP vs learning alg performance %%
strat_delta = BCRP_perf(end) - S(end);
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'expert_stats_TC','.xlsx'),stats,1) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 15. Transaction costs analysis %%
%mean TC per day 
TC_mean = mean(PnL_TC(:,1))*85*10000;
n_mean = mean(PnL_TC(:,2),'omitnan');
vol_mean = mean(PnL_TC(:,3));
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'TC_stats','.xlsx'),TC_mean,1) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'TC_stats','.xlsx'),n_mean,2) 
xlswrite(strcat(userpathstr,'\Masters\Plots\DailyDatawithTC\',savefolder,'\',string(noofportfoliostocks),'_strats',string(nostrats),'_',string(Datetimevec(1)),'-',string(Datetimevec(end)),'_',num2str(noofdays),'TC_stats','.xlsx'),vol_mean,3) 


%--------------------------------------------------------------------------