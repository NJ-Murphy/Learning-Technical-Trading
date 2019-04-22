%% Masters: Strategy testing and plotting for 2 stocks (Daily data)
%
% Author: N.J. Murphy
%
% Person Version:
% Masters-scripts-Learn2stocksDaily
%
% 1 Problem Specification:  This script aims to 
%
% 2 Data Specification: 6 months worth of order book data for 16 stocks listed
%   on the JSE Top 40. The specified period is from 01 January 2013 to 03 June
%   2013. 
%   In order to n and publish this script, 2 stocks were ralabelled in order
%   to loop over a 2 week period and 2 stocks from the test data was
%   used.
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
% Uses: This script will used to call the required strategies and plot the
% signal lines and strategy values and plot the associated buy/sell hold
% signals. 
 
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
%projectpath = 'AMF_FinalProject_NJM_715798';   %path for AMF machine
addpath(fullfile(userpathstr,projectpath,'Functions'));
addpath(fullfile(userpathstr,projectpath,'Functions','Strategies'));
addpath(fullfile(userpathstr,projectpath,'Functions','Indicators'));
addpath(fullfile(userpathstr,projectpath,'Scripts'));
addpath(fullfile(userpathstr,projectpath,'Data\TRTH\EOD'));
figpath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Plots\Strategies\MOVAVX';

% # Filepath to input data
%[datafilepath] = 'C:\Users\Nicholas\Documents\MATLAB\Masters\Data\TRTH\EOD';

% Paths for strategies
stratpath = 'C:\Users\nicjm\Documents\MATLAB\Masters\Functions\Strategies';
stratdirec = dir(fullfile(stratpath));
stratfilenames = {stratdirec(~[stratdirec.isdir]).name}';

%%%%%%%%%%%%%%%%%%
%% 5. Load Data %%
% Load data for opening and closing prices, high and low prices, volume
% traded and dates.

% Load stocks: OHLCV and choose number of stocks and number of days on from
% the 736th day

noofdays = input('Number of days = ');
%noofdays = 100;
%noofstocks = 5;
noofstocks = input('Required number of stocks = ');
%nostrats = 2;  %specify the number of strategies
%noofstrats = input('Required number of strategies = ');

% Get the letter of the column for the 1st to noofstocks in excel 
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZAAABACADAEAFAGAHAIAJAKALAMANAOAPAQARAS';
if noofstocks > 24
  stock_letter = Alphabet([2+noofstocks,3+noofstocks]);
else
  stock_letter = Alphabet(2+noofstocks);
end
stocksforexcelimport = strcat('C736:',stock_letter,num2str(736+noofdays));

% Extract the OHLCV stock data for noofstocks for the required period 
P_c = xlsread('JSEClosing',stocksforexcelimport);
P_o = xlsread('JSEOpening',stocksforexcelimport);
P_h = xlsread('JSEHigh',stocksforexcelimport);
P_l = xlsread('JSELow',stocksforexcelimport);
P_v = xlsread('JSEVol',stocksforexcelimport);

% Get the names of the stocks
[~,stocknames,~] = xlsread('JSEClosing',strcat('C2:',strcat(stock_letter,'2')));

% Get RIC for stocks
for i=1:noofstocks
    RIC = stocknames{1,i}(1:end-3);
end

% Load dates
[~,~,Dates] = xlsread('JSEClosing',strcat('A736:','A',num2str(736+noofdays)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Process and Prepare the Data for learning %%
% Here, we prepare the data for input into the learning algorithm

% Remove any stocks for which there are more than 80% of the data point
% missing
[r,c] = size(P_v); 

for i = 1:c
  if sum(isnan(P_c(:,i)))/r > 0.4  %If more than 40% of the column is NAN then remove this column
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

% Aggregate opening and closing prices, high and low prices, volume
% traded and dates into one single array (Data = [Stocks,Prices,Dates]).
Data = nan(noofstocks,5,r);
for i = 1:noofstocks
  Data(i,:,:) = [P_c(:,i)'; P_o(:,i)'; P_h(:,i)'; P_l(:,i)'; P_v(:,i)'];
end

clearvars P_c P_o P_h P_l P_v y z i c r stratdirec 

% Remove NaNs from data (weekends and public holidays)
Dates = Dates(1:size(Data,3)); %remove dates which 'Data' did not import due to empty cells
Dates(any(any(isnan(Data),1),2)) = [];
Data(:,:,any(any(isnan(Data),1),2)) = [];

% Convert data to log data
Data = log(Data);

%number of actual trading days
nooftradingdays = size(Data,3); 

% Specify k and ell values to be used as parameters and the number of
% strategies
ell = 10:14;%:30; %:12; %28;
k = 13:26;%:44; %:15; %31;

%strats = cell(3,3);
strats = stratfilenames(:,1);

% Determine the names of input parameters to see how many
% free parameters are required
for i = 1:size(strats,1)
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
[temp, order] = sort(cell2mat(strats(:,2)),'descend');
strats = strats(order,:);
store = nan(size(strats,1),length(ell),length(k));
tmin0 = max(k) + max(ell)+1;  
expind = 1:sum(cell2mat(strats(:,2)))*(length(k)*(min(k)-min(ell))+(length(k)-1)*((length(k))/2)-((max(k)-max(ell))-1)*(max(k)-max(ell))/2); %indices for experts with 2 params
expind0 = 0;  %initialise

for w0 = 1:size(strats,1)  % strategies                
    %----> Find the name of the w-th strategy
    if strcmp('RSIstrat.m', strats{w0,1})== 1 %strcmp('MACDstrat002.m', strats{w0,1}) == 1||  ||...
            % strcmp('MovAveXover.m', strats{w0,1})==1
       strategyfn = str2func(strats{w0,1}(1:end-2));  %convert the w-th strategy into a function
    else
       continue
    end
    
    for ell0 = 1:length(ell)  %parameter 1 -> n/n1
        ell1 = ell(ell0);   % ----> p.ell will be a vector of chosen lookback periods to try for n/n1 for each strategy
        for k0 = 1:length(k)   %parameter 2 -> n2/lambda    
            %----> Condition to break out of loop if only one param
            % is required for current strategy
            if cell2mat(strats(w0,2)) == 0
                break
            end                       
            k1 = k(k0); % ----> p.k will be a vector of chosen lookback periods to try for n2 for each strategy

            %----> Continue loop if the k parameter is < or =
            %to ell since ell==n/n1 cannot be > or = to k
            if k1 > ell1

               %index of the w-th expert for the current param combo
               expind0 = expind0+1;
               Expind = expind(expind0);
               
               % Store the data output for the current strategy and current set of
               % parameters
               store = strat_visualise(Data,strategyfn,strats{w0,3},ell1,k1);
               
               sig = cell(1,size(Dates,1));
               % Get signals for current combo of params 
               for t = tmin0:nooftradingdays    %time loop
                   xt = Data(:,:,1:t);
                   sig{t} = SignalGen(xt,strategyfn,strats{w0,3},ell1,k1);
                   %signals = cell2mat(sig);
               end
               
               % fill empty cells
               sig(cellfun('isempty',sig)) = {NaN};
               sigs = cell2mat(sig');
               Datetimevec = datetime(Dates);
               buysigs = Datetimevec(sigs==1);
               sellsigs = Datetimevec(sigs==-1);
               holdsigs = Datetimevec(sigs==0);
               
               % Create a weights vector to show which stock is being held
               % until a change in signal
               indnan = find(isnan(sigs),1,'last'); %find index of last NaN value
               wsigs = sigs; %weighted signals
               for i = 1:length(sigs)
                  if sigs(i)== 0 && isempty(find(sigs(indnan+1:i),1,'last'))==1
                     wsigs(i) = NaN;
                  elseif sigs(i)== 0 && isempty(find(sigs(indnan+1:i),1,'last'))~=1
                     wsigs(i) = sigs(find(sigs(indnan+1:i),1,'last')+indnan);
                  else 
                     continue
                  end
               end
               wsigs(1) = 0;
               
               % Get the required stategy name
               stratname = strats{w0,1};
               n = sub2ind([size(strats,1),length(ell),length(k)],w0,ell0,k0);
               
               switch stratname
 
                  case {'MACDstrat002.m'}
                  %figure(Expind)
                  figure;
                   %subplot(3,1,1)
              
                    plot(Datetimevec,store{1,2},Datetimevec,store{1,3});
                    datetick('x','dd-mmm-yyyy','keepticks','keeplimits')
                    xlim([Datetimevec(1) Datetimevec(end)])
                    %xticks([Datetimevec(1) Datetimevec(40) Datetimevec(70) Datetimevec(100) Datetimevec(130) Datetimevec(160) Datetimevec(end)])
                    %keeps lines starting and ending on axes
                    axis tight
                     hold on
                     set(gca,'XTickLabelRotation',50)
                     vline(buysigs(:),'g','')
                     vline(sellsigs(:),'r','')
                     %plot(x(sigs==1),store{1,2}(1,sigs==1),x(sigs==1),store{1,3}(1,sigs==1),'b.','MarkerSize',12,'MarkerFaceColor','r') % marking the 10th data point of x and y
                     
                     title1 = strcat(RIC,{': '},{'MACD and MACD Signal '},{'($\ell$='},num2str(ell1),{', $k$='},num2str(k1),{')'});
                     title(title1,'interpreter','latex')  
                     legend('MACD Signal','MACD','Location','Best')
                     xlabel('Date')
                     hold off

                     
                    figure;
                   %subplot(3,1,2) 
                    stem(Datetimevec,wsigs);
                    datetick('x','dd-mmm-yyyy','keepticks')
                    %xticks(Datetimevec(1):Datetimevec(end))
                    %xlim([Datetimevec(1) Datetimevec(end)])
                      
                    %keeps lines starting and ending on axes
                    %axis tight
 
                     hold on
                     set(gca,'XTickLabelRotation',50)
                     title1 = strcat({'Buy and Sell Weights'});
                     title(title1,'interpreter','latex')  
                     %legend('MACD Signal','MACD','Location','Best')
                     xlabel('Date')
                     hold off  
                     
                     figure;
                    %subplot(3,1,3)
                     candle(reshape(Data(:,3,:),size(Data(:,2,:),1),size(Data(:,2,:),3))',...
                         reshape(Data(:,4,:),size(Data(:,3,:),1),size(Data(:,3,:),3))',...
                         reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3))',...
                         reshape(Data(:,2,:),size(Data(:,1,:),1),size(Data(:,1,:),3))','k',Datetimevec);
                     datetick('x','dd-mmm-yyyy','keepticks')
                     ch = get(gca,'children');
                     set(ch(1),'FaceColor','r')
                     set(ch(2),'FaceColor','g')
                     
                     hold on 
                     set(gca,'XTickLabelRotation',50)                    
                     title1 = strcat({'Log of Closing Price'});
                     title(title1,'interpreter','latex')  
                     %legend('SMA Short','SMA Long','Location','Best')
                     xlabel('Date')
                     ylabel('Price ($log(P^{C}$))','interpreter','latex')
                     hold off
                    
                      
                   case {'MovAveXover.m'}   
                     %figure(Expind)
                     figure;
                     %subplot(2,1,1)
                     candle(reshape(Data(:,3,:),size(Data(:,2,:),1),size(Data(:,2,:),3))',...
                         reshape(Data(:,4,:),size(Data(:,3,:),1),size(Data(:,3,:),3))',...
                         reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3))',...
                         reshape(Data(:,2,:),size(Data(:,1,:),1),size(Data(:,1,:),3))','k',Datetimevec,'dd-mmm-yyyy');
                     ch = get(gca,'children');
                     set(ch(1),'FaceColor','r')
                     set(ch(2),'FaceColor','g')
                     datetick('x','dd-mmm-yyyy','keepticks')
                     %datetick('x','dd-mmm-yyyy','keepticks')
                     %axis tight

                     hold on 
                     h = plot(datenum(Datetimevec),store{1,2},datenum(Datetimevec),store{1,3});

                     vline(datenum(buysigs(:)),'g','')
                     vline(datenum(sellsigs(:)),'r','')
                     %plot(x(sigs==1),store{1,2}(1,sigs==1),x(sigs==1),store{1,3}(1,sigs==1),'b.','MarkerSize',12,'MarkerFaceColor','r') % marking the 10th data point of x and y
                     set(gca,'XTickLabelRotation',50)
                     
                     title1 = strcat(RIC,{': '},{'Moving Avg Crossover ('},{'$\ell$='},num2str(ell1),{', $k$='},num2str(k1),{')'});
                     title(title1,'interpreter','latex')  
                     legend(h,'SMA Short','SMA Long','Location','Best')
                     xlabel('Date')
                     ylabel('Price ($log(P^{C}$))','interpreter','latex')
                     hold off

                     %figure;
                     %subplot(2,1,2)

                    % hold on 
                    % set(gca,'XTickLabelRotation',50)                    
                    % title1 = strcat({'Log of Closing Price'});
                    % title(title1,'interpreter','latex')  
                     %legend('SMA Short','SMA Long','Location','Best')
                    % xlabel('Date')
                    % ylabel('Price ($log(P^{C}$))','interpreter','latex')
                    % hold off
                     
               end
               clearvars store
            else 
                continue
            end 
        end % k

        %---> Compute the agents performance and find the optimal
        % parameters if the strategy only requires one
        % parameter
        if cell2mat(strats(w0,2)) == 0
            expind0 = expind0+1;
            % Store the data output for the current strategy and current set of
               % parameters
               store = strat_visualise(Data,strategyfn,strats{w0,3},ell1,{});
               
               sig = cell(1,size(Dates,1));
               % Get signals for current combo of params 
               for t = tmin0:nooftradingdays    %time loop
                   xt = Data(:,:,1:t);
                   sig{t} = SignalGen(xt,strategyfn,strats{w0,3},ell1,{});
                   %signals = cell2mat(sig);
               end
               
               % fill empty cells
               sig(cellfun('isempty',sig)) = {NaN};
               sigs = cell2mat(sig');
               Datetimevec = datetime(Dates);
               buysigs = Datetimevec(sigs==1);
               sellsigs = Datetimevec(sigs==-1);
               holdsigs = Datetimevec(sigs==0);
               
               %check if there are cases where there are no buy/sell
               %signals
               if isempty(sellsigs)==1
                   sellsigs = NaT;
               elseif isempty(buysigs)==1
                   buysigs = NaT;
               end
               
               % Get the required stategy name
               stratname = strats{w0,1};
           
            switch stratname
              case {'RSIstrat.m'}
               figure;
                  % figure(expind0)
                 %subplot(2,1,1)
                 plot(Datetimevec,store{1,2});
                 %axis tight

                 hold on  
                 %hax1 = get(h1, 'Parent'); % axes handle
                 %set(gca,'XTickLabelRotation',50)
                 %set(hax1, 'XTick', dateV(1:15:end)); 
                 datetick('x','dd-mmm-yyyy','keepticks','keeplimits')
                 hline([30 70],{'k','k'},{'Oversold (30)','Overbought (70)'})
                 vline(buysigs(:),'g','')
                 vline(sellsigs(:),'r','')
                 %plot(x(sigs==1),store{1,2}(1,sigs==1),x(sigs==1),store{1,3}(1,sigs==1),'b.','MarkerSize',12,'MarkerFaceColor','r') % marking the 10th data point of x and y
                 set(gca,'XTickLabelRotation',50)

                 title1 = strcat(RIC,{': '},{'Relative Strength Index'},{' ($\ell$='},num2str(ell1),{')'});
                 title(title1,'interpreter','latex')  
                 legend('RSI','Location','Best')
                 xlabel('Date')
                 ylabel('RSI')
                 hold off

                
                 figure;
                 %subplot(2,1,2)
                 candle(reshape(Data(:,3,:),size(Data(:,2,:),1),size(Data(:,2,:),3))',...
                     reshape(Data(:,4,:),size(Data(:,3,:),1),size(Data(:,3,:),3))',...
                     reshape(Data(:,1,:),size(Data(:,1,:),1),size(Data(:,1,:),3))',...
                     reshape(Data(:,2,:),size(Data(:,1,:),1),size(Data(:,1,:),3))','k',Datetimevec);
                 ch = get(gca,'children');
                 set(ch(1),'FaceColor','r')
                 set(ch(2),'FaceColor','g')
                 datetick('x','dd-mmm-yyyy','keepticks')
                 hold on 
                 set(gca,'XTickLabelRotation',50)                    
                 title1 = strcat({'Log of Closing Price'});
                 title(title1,'interpreter','latex')  
                 %legend('SMA Short','SMA Long','Location','Best')
                 xlabel('Date')
                 ylabel('Price ($log(P^{C}$))','interpreter','latex')
                 hold off
                  
                                 %%%%%%%%%%%%%%%
                % save figure %
                %%%%%%%%%%%%%%%
                savetoname = char(strcat(figpath,'\',num2str(noofdays),RIC,'logprices'));
                fig = gcf;
                fig.PaperUnits = 'centimeters';
                %fig.PaperPosition = [0 0 19 12]; %proposal size
                fig.PaperPosition = [0 0 12 12];
                print(savetoname,'-dpng','-r300')
                %case {'Bollstrat'}   

           end
        end
    end % ell
end % w   

%--------------------------------------------------------------------------   