classdef learning_IntradayDaily
    % The class implements online technical trading and
    % learning over time for M objects (stocks) and P features (OHLCV) as specified in
    % an MxPxT data matrix X. The algorithm computes the trading signals and associated
    % for each expert at each time period. An expert is characterised by a
    % given combination of the 4 free parameters in the algorithm, namely,
    % w, c, k and ell corresponding to the strategy, cluster, long-term look-back
    % parameter and short-term look-back parameter. Using the controls H(W,K,L,CI;T)
    % the expert wealth SH(K,L,CI;T) is computed. The controls are then aggregate 
    % using machine learning to provide the overall portfolio controls b and 
    % realised cumulative overall portfolio wealth S.
    %
    % See Also: PATTERN/CONTROLS, PATTERN/OFFLINE, PATTERN/LEARN, PATTERN/ONLINE
    %
    %
    % A. OHLC patterns
    % B. Fundamental model patterns
    %
    % 1. *Data* 
    %   1.1 M  entities (objects - stocks)
    %   1.2 P  features (OHLCV)
    %   1.3 T  observations (date-times)
    % 2. *Parametrs*
    %   2.1 k vector
    %   2.2 ell vector
    %   2.3 set of W strategies
    %   2.4 clusters vector
    % 3. *Experts*
    %   3.1. Active
    % 4. *Algorithm*
    %   4.1. Online/Offline
    %   4.2. Interface experts (h,SH)
    % 5. *Machine Learning* (UNIV + active)
    %   5.1. Weighted Arithmetic Average over all experts to create optimal predictors
    %       5.1.1. Performance weighted averaging (using arithmetic averaging)
    %           B(T) = SUM(K,L) (SH(T-1|K,L) H(T|K,L)) / SUM(K,L) SH(T-1|K,L)
    %           this is probability wieghted where the probability is
    %           proportional to the returns.
    % 6. *Sectors and States* (clusters)
    %   6.1. CxP for C clusters and M objects
    
    % Authors: N. Murphy, T. Gebbie
    
    %% public properties
    properties
        b = []; % aggregated controls    (T,M) Time x Objects  ---- best combined strategies 
        S = []; % aggregated performance (T,1) Time x 1  -----  weights 
        h = []; % agent controls      (N,M,T)  Agents x Objects x Time  ---- 
        H = []; %daily agent controls
        SH_fused = []; % agents performance (T,N)    Time x Agents --- Trading Strategy- time series/ strategies/ rules   
        SHave_fused = []; %mean wealth for each strategy 
        PnL = []; %profit and losses overall strategy
        PnL_experts = []; %profit and losses experts
        %expertparams = []; %stores the 4 free parameters used in the algorithm
        %stockstraded = [];  %matrix to keep track of 'noofportfolio' most liquid stocks every periodinterval periods
        S_TC = [];
        PnL_TC = []; 
        TC_vec = [];
    end
    
    %% private properties
    properties(Access = private)
        %x = []; % data as price relatives (M,F,T) Objects x Features x Time
        x_daily = [];
        x_intraday = [];
        uniquedaytr = [];
        k = []; % vector of values for parameter 2- n2
        ell = []; % vector of values for parameter 1- n/n1
        p = []; % current partition
        dimIntraday = []; %Intraday data dimensions
        dimdaily = []; %daily data dimensions
        ntype = 'active'; % agent normalisation (estimation)
        lntype = 'active'; % agent relative (learning) normalisation
        ltype = 'univlearning';  
        lparam = []; % learning parameters
        ptype = 'trivial'; %partition type
        weighttype = ''; %signal to weight transformation type (volatility or inverse volatility)
        qH = []; % agents controls (T,N) Time x Agents
        strats = {}; % names of the strategies
        STEFI = []; %STEFI index 
        sectors = {}; %sector clusters
        tmin0 = []; %number of periods before online learning begins
        nostrats = [];  % number of strategies
        indstorel = [];
        indstorekl = [];
        periodinterval = []; %number of days after liquidity of stocks is checked
        stocknames = []; %names of stocks needed to find clusters
        noofportfoliostocks = []; %number of stocks to be used in the portfolio construction throughout time
        liquidstocks = []; %indices of 'noofportfoliostocks' most liquid stocks at each periodinterval
        hlt = [];  %contols for strategies with 1 parameter
        hlkt = [];  %contols for strategies with 2 parameters
        %fused parameters (daily and intraday)
        noofexperts = [];
        hlkt_daily = [];
        hlt_daily = [];
    end
    
    %% methods
    methods
        %% constructor
        function p = learning_IntradayDaily(varargin)
            % P = LEARNING/LEARNING constructor
            %
            % P = LEARNING(X,K,ELL,CI,NTYPE,LNTYPE,LTYPE) For data X (price relatives)
            %   and is MxNxP dimensional data for M objects, N features and
            %   P date-times. dim(X) is dim(Price)-1. A typical object is
            %   a stocks, e.g. AGL, a typical feature is a factors such as
            %   OPEN, HIGH, LOW, CLOSE (OHLC), and date-times are the time
            %   stamps for the data.
            %
            % See Also LEARNING/DISPLAY, LEARNING/PARTITION
            
            if nargin==0
                % null constructor
            else
                %% Compute price relatives as returns
                % Data
                x_intraday = varargin{1}; %this corresponds to the input Data
                x_daily = varargin{2};
                
                % size of the data (m=stocks, n=features, p=dates/times)
                dimIntraday = size(x_intraday);
                dimdaily = size(x_daily);
                
                params2 = varargin{3}; %extract the parameter values for k
                params1 = varargin{4}; %extract the parameter values for ell
                
                %---> specify number of strategies to be used
                strats = varargin{5}; 
                
                %---> input STEFI
                STEFI = varargin{6};
                
                periodinterval = varargin{7};
                
                % ---> stock names
                stocknames = varargin{8};
                
                %---> sector clusters
                sectors = varargin{9};
                
                %---> stocks for portfolio contruction
                noofportfoliostocks = varargin{10};
                
                liquidstocks = varargin{11};
                
                weighttype = varargin{12};
                
                uniquedaytr = varargin{13};
                
                %---> set defaults for optional inputs
                optargs = {x_intraday x_daily params2 params1 strats STEFI periodinterval stocknames...
                    sectors noofportfoliostocks liquidstocks weighttype uniquedaytr};

                % Place optional args in memorable variable names
                [p.x_intraday, p.x_daily, p.k, p.ell, p.strats, p.STEFI, p.periodinterval,...
                    p.stocknames, p.sectors, p.noofportfoliostocks,...
                    p.liquidstocks, p.weighttype, p.uniquedaytr] = optargs{:};
                
                % update the size property
                p.dimIntraday = dimIntraday;
                p.dimdaily = dimdaily;
                
                %---> Number of sectors
                C = length(p.sectors);
                
                %---> control parameters
                L = size(p.ell,2); % number of param 1
                
                % number of experts 1 paramter and 2 parameters
                noofexperts1param = C*L*sum(cell2mat(strats(:,2))==0);
                noofexperts2params = C*sum(cell2mat(p.strats(:,2)))*sum(sum(p.k>max(p.ell)):sum(p.k>min(p.ell)));
                p.noofexperts = noofexperts2params + noofexperts1param;
                
                % initialise state variables and pre-allocate memory
                % Intraday from initial to terminal time T:
                p.b = nan(p.dimIntraday(3)+p.dimdaily(3),noofportfoliostocks+1); % Time x Objects (assets + rf asset)   
                p.h = nan(p.noofexperts,noofportfoliostocks+1,p.dimIntraday(3)+2); % Agents x Objects(stocks + rf asset) x Time    
                p.qH = nan(p.dimIntraday(3)+p.dimdaily(3)+1,p.noofexperts); % Time x Experts
                p.S = ones(p.dimIntraday(3)+p.dimdaily(3),1); % Time x 1
                %p.SH = ones(p.dimIntraday(3)+2,noofexperts); %  Time x Experts  (length(p.sectors) for the stock clusters
                %p.SHave = ones(p.dimIntraday(3)+2,size(strats,1)); %Time x strategies
                p.hlkt = zeros(noofexperts2params,p.noofportfoliostocks+1,p.dimIntraday(3)+2); %matrix for strats with 2 params
                p.hlt = zeros(noofexperts1param,p.noofportfoliostocks+1,p.dimIntraday(3)+2); %matrix with strats for 1 param
                p.PnL = zeros(p.dimIntraday(3)+p.dimdaily(3),1); % Time x 1
                p.S_TC = ones(p.dimIntraday(3)+p.dimdaily(3),1); % Time x 1
                p.PnL_TC = zeros(p.dimIntraday(3)+p.dimdaily(3),1); % Time x 1
                p.PnL_experts = zeros(p.dimIntraday(3)+p.dimdaily(3),p.noofexperts); %  Time x Experts
                p.TC_vec = zeros(p.dimIntraday(3)+p.dimdaily(3),3); % Time x 3

                %Daily:
               % p.B = nan(p.dimIntraday(3)+p.dimdaily(3)+1,noofportfoliostocks+1); % Time x Objects (assets + rf asset)   
                p.H = nan(p.noofexperts,noofportfoliostocks+1,p.dimIntraday(3)+p.dimdaily(3)+1); % Agents x Objects(stocks + rf asset) x Time    
                %p.qH_fused = nan(p.dimIntraday(3)+p.dimdaily(3)+1,noofexperts); % Time x Experts
                %p.S_fused = ones(p.dimIntraday(3)+p.dimdaily(3)+1,1); % Time x 1
                p.SH_fused = ones(p.dimIntraday(3)+p.dimdaily(3),p.noofexperts); %  Time x Experts  (length(p.sectors) for the stock clusters
                p.SHave_fused = ones(p.dimIntraday(3)+p.dimdaily(3),size(strats,1)); %Time x strategies
                
                p.hlkt_daily = zeros(noofexperts2params,p.noofportfoliostocks+1,p.dimdaily(3)+1); %matrix for strats with 2 params
                p.hlt_daily = zeros(noofexperts1param,p.noofportfoliostocks+1,p.dimdaily(3)+1); %matrix with strats for 1 param
              
                %% partition the data
                p = partition_Intraday(p);
            end  
        end
        
        %% offline (to force offline estimation)
        function p = offline_Intraday(p)
            % LEARNING/OFFLINE Offline estimation
            %
            % P = OFFLINE(P) to estimate (H,SH) for T for parameters k,
            % ell and c
            %
            % See Also: LEARNING/ONLINE
            
            %---> Specify the number of periods for the offline segement of
            % learning
            p.tmin0 = 2;
            T = p.dimdaily(3); %total number of time periods in the data
            
            if T < p.tmin0
                error('pattern:offline','Not enough Data L*K>T');
            end
            
            % compute returns for periods 1 to tmin0 (before all trading
            % starts)
            % Extract closing prices to compute returns
            Pc_daily = reshape(p.x_daily(:,1,1:p.tmin0-1),size(p.x_daily(:,1,1:p.tmin0-1),1),size(p.x_daily(:,1,1:p.tmin0-1),3)); %(stocks x time)

            % Method 1: Compute arithmetic returns for closing s and returns
            %FIX ME!!
            %ret = [zeros(1,size(Pc,1)+1);[diff(Pc')./Pc(:,1:end-1)',diff(p.STEFI(1:t,1))./p.STEFI(1:(t-1),1)]]; %add row of zeros at the start of vector to account for 1st period ret
            
            % Method 2: Compute log returns
            %FIX ME!!
            %ret = [zeros(1,size(Pc,1)+1);[diff(log(Pc),1,2)',diff(log(p.STEFI(1:t,1)))]]; %add row of zeros at the start of vector to account for 1st period ret
            %ret = ret(:,[p.liquidstocks(:,1:p.noofportfoliostocks) end]);
            
            % Method 3: Compute price relatives
            ret_all_daily = [ones(1,size(Pc_daily,1)+1);[transpose(Pc_daily(:,2:end)./Pc_daily(:,1:end-1)),p.STEFI(2:p.tmin0-1,1)./p.STEFI(1:p.tmin0-1-1,1)]]; %only compute returns up to befre tmin0 as tmin0
            
            % set controls for initial time period (t=1) as no trading takes place throughout this day
            p.h(:,:,1:p.uniquedaytr(2)) = zeros;%(size(p.h(:,:,1)));
            p.b(1:p.uniquedaytr(2),:) = zeros;
            
            % implement the sequential learning for next time t
            for t = p.tmin0:T   %day loop
               %initialise intraday returns for t-th day 
               ret_all_intra = ones(1,size(Pc_daily,1)+1);
               
               p.SH_fused(p.uniquedaytr(t)+t-1) = ones;
               %p.SH_fused(p.uniquedaytr(t)+t-1) = p.SH_fused(p.uniquedaytr(t)+t-2);
               
               %since STEFI is a daily rate, keep same STEFI throughout
               %each day
               STEFI_intra = repmat(p.STEFI(t,1),p.uniquedaytr(t+1) - p.uniquedaytr(t),1);
                
               for t_intra = (p.uniquedaytr(t)+1):p.uniquedaytr(t+1)-1 %intraday loop: +2 and not -1 from p.uniquedaytr(t+1) since we want the first entry of each day in h, S, etc to be the daily trade
                   Pc_intra = reshape(p.x_intraday(:,1,t_intra-1:t_intra),size(p.x_intraday(:,1,t_intra-1:t_intra),1),size(p.x_intraday(:,1,t_intra-1:t_intra),3)); %closing prices time t-1 to t  
                   ret_all_intra = [ret_all_intra; [transpose(Pc_intra(:,end)./Pc_intra(:,end-1)),STEFI_intra(t_intra-p.uniquedaytr(t)+1,1)./STEFI_intra(t_intra-p.uniquedaytr(t),1)]];
                   p = online_Intraday(p,t_intra,t,ret_all_intra);
               end
               %Daily learning 
               Pc_daily = reshape(p.x_daily(:,1,t-1:t),size(p.x_daily(:,1,t-1:t),1),size(p.x_daily(:,1,t-1:t),3)); %closing prices time t-1 to t                
               %ret_all = [ret_all; [transpose(Pc(:,end)./Pc(:,end-1)),p.STEFI(t,1)./p.STEFI(t-1,1)]];
               ret_all_daily = [ret_all_daily; [transpose(Pc_daily(:,end)./Pc_daily(:,end-1)),p.STEFI(t,1)./p.STEFI(t-1,1)]];
               p = online_Interday(p,t,p.uniquedaytr(t+1),ret_all_daily); 
            end      
        end
        
        %% online (to force online estimation)
        function p = online_Intraday(varargin)
            % LEARNING/ONLINE Offline estimation
            %
            % P = ONLINE(P) to estimate (H,SH,B,S) at T for the range of
            %   K and L over the specified clusters CI. This requires the
            %   online structure for learning to have been initialised. T
            %   is taken to be the last time in the object.
            %
            % P = ONLINE(P,t) to estimate online values at time t using the
            %   data from times 1 to t where t = {1:T}.
            %
            % See Also: LEARNING/LEARN, LEARNING/CONTROLS
            
            %% input parameters

            %---> now put these defaults into the valuesToUse cell array,
            optargs(1:nargin) = varargin; %varargin = (p,t)
            %---> Place optional args in memorable variable names
            [p, t_intra, t, ret_all_intra] = optargs{:};
            
            %% control parameters
            L = sum(double(p.ell <= t_intra-p.uniquedaytr(t))); % number of param 1
            W = size(p.strats,1); % number of strategies
            K = sum(double(p.k <= t_intra-p.uniquedaytr(t))); % number of params 2
            
            %% Model Identification Loop
            % Extract data from time t=1 to current time
            xt = p.x_intraday(:,:,p.uniquedaytr(t):t_intra);
            
            % only get returns for the liquid stocks we are considernt for
            % current 'periodinterval' period
            ret = ret_all_intra(:,[p.liquidstocks(1:p.noofportfoliostocks) end]); %(time x stocks) 
              
            %check if sufficient data to run at least one parameter loop
            if all(p.ell >= t_intra-p.uniquedaytr(t))
                if t_intra == (p.uniquedaytr(t)+1) %case for beginning of the day when we have to skip first time period to compute a return
                    t_initial = t_intra - 1;
                    p.h(:,:,t_initial+1) = zeros; %fill nans in with zeros as no experts trade in this case
                    
                    %update intraday-daily wealth matrix
                    dSH = ones(size(p.SH_fused(t_intra,:)));
                    p.SH_fused(t_initial+t-1,:) = ones;%p.SH_fused(t_initial+t-2,:).*dSH;
                    ret_initial = ones(size(ret)); 
                    p = learn_Intraday(p,t_initial+t-2,ret_initial,p.h(:,:,t_initial+1),t_intra,t);
                end
                p.h(:,:,t_intra+1) = zeros; %fill nans in with zeros as no experts trade in this case
                dSH = ones(size(p.SH_fused(t_intra,:)));
                %p.SH(t_intra+1,:) =  p.SH(t_intra,:).*dSH;
                %update intraday-daily wealth matrix
                p.SH_fused(t_intra+t-1,:) = p.SH_fused(t_intra+t-2,:).*dSH;
                p = learn_Intraday(p,t_intra+t-2,ret,p.h(:,:,t_intra+1),t_intra,t);
                return
            end
                    
            % Initialise the signal matrices for the t-th stime step
            C = length(p.sectors); %number of clusters
            
            %initialise indices for experts
            hlktind1 = 0;  
            hltind = 0;
            expertparam2 = zeros(p.noofexperts,4);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% ---> Going to need to edit this section to learn
            %%% strategies. Instead of using the matching algorithm to
            %%% create the agents
            
            % strategy loop
            for c = 1:C % clusters (sectors) of stocks 
                for w0 = 1:W % strategies               
                    %----> Find the name of the w-th strategy
                    strategy = p.strats{w0,1}(1:end-2);
                    strategyfn = str2func(strategy);  %convert the w-th strategy into a function

                    for ell0 = 1:L  %parameter 1 -> n/n1
                        ell1 = p.ell(ell0);   % ----> p.ell will be a vector of chosen lookback periods to try for n/n1 for each strategy
                        %check if theres sufficient data for lookback param
                        %ell
                        if ell1 >= t_intra-p.uniquedaytr(t)
                            break
                        end
                        for k0 = 1:K   %parameter 2 -> n2/lambda  
                            k1 = p.k(k0); % ----> p.k will be a vector of chosen lookback periods to try for n2 for each strategy

                            if k1 > t_intra-p.uniquedaytr(t)
                                break
                            end
                            %----> Condition to break out of loop if only one param
                            % is required for current strategy
                            if cell2mat(p.strats(w0,2)) == 0
                                break
                            end                       

                            %----> Continue loop if the k parameter is < or =
                            %to ell since ell==n/n1 cannot be > or = to k
                            if k1 > ell1
                               %Expert index KLrow(w0,k0,ell0)
                               hlktind1 = hlktind1+1;
                               Expind = hlktind1;
                               %Expind = hlktind(hlktind1); 
                               p.indstorekl(Expind) = w0;
                               %Compute the signals of the w-th expert for the
                               %current combination of k and ell
                               sig = controls_Intraday(xt(logical(p.sectors{c}),:,:),strategyfn,p.strats{w0,3},p.hlkt(:,:,t_intra),Expind,p.weighttype,ell1,k1); %(logical(p.sectors{c}),:,:)
                               p.hlkt(Expind,p.sectors{c}>0,t_intra+1) = sig(1:end-1);  %Experts x stocks in cluster c
                               p.hlkt(Expind,sum(p.sectors{c}>0)+1,t_intra+1) = sig(end); %risk free weight
                               expertparam2(Expind,:) = [c,w0,ell1,k1];
                            else 
                                continue
                            end 
                        end % k

                        %---> Compute the agents performance and find the optimal
                        % parameters if the strategy only requires one
                        % parameter
                        if cell2mat(p.strats(w0,2)) == 0
                           %this will contain the same information as the if
                           %statement in the k loop above
                           hltind = hltind + 1;
                           Expind = hltind;
                           
                           p.indstorel(Expind) = w0; %store the strategy indices

                           % check if its anticor, bcrp or antbcrp:
                           if func2str(strategyfn) == "OnlineZAnticor"
                              sig = OnlineZAnticor(ret(:,logical([p.sectors{c},1])),size(ret,1), ell1, p.hlt(Expind,logical([p.sectors{c},1]),t_intra));
                              p.hlt(Expind,logical([p.sectors{c}>0,1]),t_intra+1) = sig;  %Experts x stocks in cluster c
                           elseif func2str(strategyfn) == "OnlineZBCRP"
                              sig = OnlineZBCRP(ret(:,logical([p.sectors{c},1])), ell1);
                              p.hlt(Expind,logical([p.sectors{c}>0,1]),t_intra+1) = sig;  %Experts x stocks in cluster c
                           elseif func2str(strategyfn) == "OnlineAntiZBCRP"
                              sig = OnlineAntiZBCRP(ret(:,logical([p.sectors{c},1])), ell1);
                              p.hlt(Expind,logical([p.sectors{c}>0,1]),t_intra+1) = sig;  %Experts x stocks in cluster c
                           else 
                              sig = controls_Intraday(xt(logical(p.sectors{c}),:,:),strategyfn,p.strats{w0,3},p.hlt(:,:,t_intra),Expind,p.weighttype,ell1);  %(logical(p.sectors{c}),:,:)
                              %aggregate current experts controsl into H
                              p.hlt(Expind,p.sectors{c}>0,t_intra+1) = sig(1:end-1);  %Experts x stocks in cluster c
                              p.hlt(Expind,sum(p.sectors{c})+1,t_intra+1) = sig(end); %risk free weight
                           end
                        end
                    end % ell                
                end % w          
            end % c
            
            % Compute the update performance for the prior agent step
            % --------------->>> This is step 2 of the online algo (pg. 4 of
            %Gebbie et al. (2016))
            if any(isnan(p.h(:,:,t_intra))) 
               dSH = ones(size(p.SH(t_intra,:)));  %first iteration
            else
               dSH = (p.h(:,:,t_intra)*(ret(end,:)'-1) + 1)';  %ret = price relatives
               %dSH = (p.h(:,:,t)*ret_clean(end,:)' + 1)';  %ret = log returns
            end
            
            %check for case where no experts using 2 params traded during
            %current period
            if sum(all(expertparam2== 0))>0
                %expertparam2 = zeros(hlktind1,4);
                hlkt1 = zeros(size(p.h(:,:,t_intra+1),1)-size(p.hlt(:,:,t_intra+1),1),size(p.hlt(:,:,t_intra+1),2));
            else
               %expertparam2 = expertparam2(1:hlktind1,:);
               hlkt1 = p.hlkt(1:hlktind1,:,t_intra+1); %only pick out rows for which k loop was activated-this is included due            
               %to the fact that we used increments of 3 in k and ell in
               %script file
               if size([hlkt1;p.hlt(:,:,t_intra+1)],1) ~= size(p.h(:,:,t_intra+1),1) %- size(p.hlt(:,:,t_intra+1),1))
                   hlkt1 = [hlkt1;zeros(size(p.h(:,:,t_intra+1),1)-size(p.hlt(:,:,t_intra+1),1)-size(hlkt1,1),size(p.hlt(:,:,t_intra+1),2))];
               end
            end
            
            % Update the agents at time t+1 and calculate their corresponding wealth
            p.h(:,:,t_intra+1) = [hlkt1;p.hlt(:,:,t_intra+1)]; %concatenate experts with 1 and 2 params resp.

            % Update the agent accumulated performance (geometric returns)
            %------->>> Here SH is the agent wealth- also need to take into
            % account borrowing costs so subtract STEFI when shorting risk free asset 
            p.SH_fused(t_intra+t-1,:) = p.SH_fused(t_intra+t-2,:).*dSH;
            
            % profits and losses of experts
            p.PnL_experts(t_intra+t-1,:) = dSH-1;
            
            %% online update the learning
            p = learn_Intraday(p,t_intra+t-2,ret,p.h(:,:,t_intra+1),t_intra,t);
        end
        
        function p = online_Interday(varargin)
            % LEARNING/ONLINE Offline estimation
            %
            % P = ONLINE(P) to estimate (H,SH,B,S) at T for the range of
            %   K and L over the specified clusters CI. This requires the
            %   online structure for learning to have been initialised. T
            %   is taken to be the last time in the object.
            %
            % P = ONLINE(P,t) to estimate online values at time t using the
            %   data from times 1 to t where t = {1:T}.
            %
            % See Also: LEARNING/ONLINE, LEARNING/CONTROLS
            
            %% input parameters

            %---> now put these defaults into the valuesToUse cell array,
            optargs(1:nargin) = varargin; 
            %---> Place optional args in memorable variable names
            [p, t, t_endofday,ret_all_daily] = optargs{:};
            
            %% control parameters
            L = sum(double(p.ell <= t)); % number of param 1
            W = size(p.strats,1); % number of strategies
            K = sum(double(p.k <= t)); % number of params 2
            
            %% Model Identification Loop
            % initial controls
            % t+1 rows as the last row is for an unrealised return
            
            % Extract data from time t=1 to current time
            xt = p.x_daily(:,:,1:t);
            
            % only get returns for the liquid stocks we are considernt for
            % current 'periodinterval' period
            ret = ret_all_daily(:,[p.liquidstocks(:,1:p.noofportfoliostocks) end]); %(time x stocks) 
            
            %check if sufficient data to run at least one parameter loop
            if all(p.ell >= t) %= is there since Bollstrat needs more than min(ell) data points
                p.H(:,:,t+1) = zeros; %fill nans in with zeros as no experts trade in this case
                dSH = ones(size(p.SH_fused(t_endofday,:)));
                %p.SH(t_endofday+1,:) =  p.SH(t_endofday,:).*dSH;
                p.SH_fused(t_endofday+t-1,:) = p.SH_fused(t_endofday+t-2,:).*dSH;
                p = learn_Intraday(p,t_endofday+t-2,ret,p.H(:,:,t+1),t_endofday,t);
                return
            end
            
            % Initialise the signal matrices for the t-th stime step
            C = length(p.sectors); %number of clusters

            %initialise indices for experts
            hlktind1 = 0; 
            hltind = 0;
            expertparam2 = zeros(p.noofexperts,4);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% ---> Going to need to edit this section to learn
            %%% strategies. Instead of using the matching algorithm to
            %%% create the agents
            
            % strategy loop
            for c = 1:C % clusters (sectors) of stocks 
                for w0 = 1:W % strategies               
                    %----> Find the name of the w-th strategy
                    strategy = p.strats{w0,1}(1:end-2);
                    strategyfn = str2func(strategy);  %convert the w-th strategy into a function

                    for ell0 = 1:L  %parameter 1 -> n/n1
                        ell1 = p.ell(ell0);   % ----> p.ell will be a vector of chosen lookback periods to try for n/n1 for each strategy
                        
                        %check if theres sufficient data for lookback param
                        %ell
                        if ell1 >= t
                            break
                        end
                        for k0 = 1:K   %parameter 2 -> n2/lambda    
                            %----> Condition to break out of loop if only one param
                            % is required for current strategy
                            if cell2mat(p.strats(w0,2)) == 0
                                break
                            end                       
                            k1 = p.k(k0); % ----> p.k will be a vector of chosen lookback periods to try for n2 for each strategy
                            
                            if k1 > t
                                break
                            end

                            %----> Continue loop if the k parameter is < or =
                            %to ell since ell==n/n1 cannot be > or = to k
                            if k1 > ell1
                               %Expert index KLrow(w0,k0,ell0)
                               hlktind1 = hlktind1+1;
                               Expind = hlktind1; 
                               p.indstorekl(Expind) = w0;
                               %Compute the signals of the w-th expert for the
                               %current combination of k and ell
                               sig = controls_Intraday(xt(logical(p.sectors{c}),:,:),strategyfn,p.strats{w0,3},p.hlkt_daily(:,:,t),Expind,p.weighttype,ell1,k1); %(logical(p.sectors{c}),:,:)
                               p.hlkt_daily(Expind,p.sectors{c}>0,t+1) = sig(1:end-1);  %Experts x stocks in cluster c
                               p.hlkt_daily(Expind,sum(p.sectors{c}>0)+1,t+1) = sig(end); %risk free weight
                               expertparam2(Expind,:) = [c,w0,ell1,k1];
                            else 
                                continue
                            end 
                        end % k

                        %---> Compute the agents performance and find the optimal
                        % parameters if the strategy only requires one
                        % parameter
                        if cell2mat(p.strats(w0,2)) == 0
                           %this will contain the same information as the if
                           %statement in the k loop above
                           hltind = hltind + 1;
                           Expind = hltind;
                           
                           p.indstorel(Expind) = w0; %store the strategy indices

                           % check if its anticor, bcrp or antbcrp:
                           if func2str(strategyfn) == "OnlineZAnticor"
                              sig = OnlineZAnticor(ret(:,logical([p.sectors{c},1])),size(ret,1), ell1, p.hlt_daily(Expind,logical([p.sectors{c},1]),t));
                              p.hlt_daily(Expind,logical([p.sectors{c}>0,1]),t+1) = sig;  %Experts x stocks in cluster c
                           elseif func2str(strategyfn) == "OnlineZBCRP"
                              sig = OnlineZBCRP(ret(:,logical([p.sectors{c},1])), ell1);
                              p.hlt_daily(Expind,logical([p.sectors{c}>0,1]),t+1) = sig;  %Experts x stocks in cluster c
                           elseif func2str(strategyfn) == "OnlineAntiZBCRP"
                              sig = OnlineAntiZBCRP(ret(:,logical([p.sectors{c},1])), ell1);
                              p.hlt_daily(Expind,logical([p.sectors{c}>0,1]),t+1) = sig;  %Experts x stocks in cluster c
                           else 
                              sig = controls_Intraday(xt(logical(p.sectors{c}),:,:),strategyfn,p.strats{w0,3},p.hlt_daily(:,:,t),Expind,p.weighttype,ell1);  %(logical(p.sectors{c}),:,:)
                              %aggregate current experts controsl into H
                              p.hlt_daily(Expind,p.sectors{c}>0,t_endofday+1) = sig(1:end-1);  %Experts x stocks in cluster c
                              p.hlt_daily(Expind,sum(p.sectors{c})+1,t_endofday+1) = sig(end); %risk free weight
                           end
                           % store the expert parameters for the current
                           % expert
                           %expertparam1(Expind,:) = [c,w0,ell1,nan];
                        end
                    end % ell                
                end % w          
            end % c
            
            % Compute the update performance for the prior agent step
            % --------------->>> This is step 2 of the online algo (pg. 4 of
            %Gebbie et al. (2016))
            if any(isnan(p.H(:,:,t))) 
               dSH = ones(size(p.SH(t,:)));  %first iteration
            else
               dSH = (p.H(:,:,t)*(ret(end,:)'-1) + 1)';  %ret = price relatives
               %dSH = (p.h(:,:,t)*ret_clean(end,:)' + 1)';  %ret = log returns
            end
      
             if sum(all(expertparam2 == 0))>0
                %expertparam2 = expertparam2(1:hlktind1,:);
                hlkt1 = zeros(size(p.H(:,:,t+1),1)-size(p.hlt_daily(:,:,t+1),1),size(p.hlt_daily(:,:,t+1),2));
            else
               %expertparam2 = expertparam2(1:hlktind1,:);
               hlkt1 = p.hlkt_daily(1:hlktind1,:,t+1); %only pick out rows for which k loop was activated-this is included due            
               %to the fact that we used increments of 3 in k and ell in
               %script file
               if size([hlkt1;p.hlt_daily(:,:,t+1)],1) ~= size(p.H(:,:,t+1),1)
                   hlkt1 = [hlkt1;zeros(size(p.H(:,:,t+1),1)-size(p.hlt_daily(:,:,t+1),1)-size(hlkt1,1),size(p.hlt_daily(:,:,t+1),2))];
               end
             end 
            
            % Update the agents at time t+1 and calculate their corresponding wealth
            p.H(:,:,t+1) = [hlkt1;p.hlt_daily(:,:,t+1)]; %concatenate experts with 1 and 2 params resp. interday controls
            
            % Update the agent accumulated performance (geometric returns)
            %------->>> Here SH is the agent wealth- also need to take into
            % account borrowing costs so subtract STEFI when shorting risk free asset 
            %p.SH(t_endofday+1,:) =  p.SH(t_endofday,:).*dSH;
            p.SH_fused(t_endofday+t-1,:) = p.SH_fused(t_endofday+t-2,:).*dSH;
            
            % profits and losses of experts
            p.PnL_experts(t_endofday+t-1,:) = dSH-1;
            
            %% online update the learning
            p = learn_Intraday(p,t_endofday+t-2,ret,p.H(:,:,t+1),t_endofday,t);
        end
        
        %% universal learning
        function p = learn_Intraday(varargin)
            % LEARNING/LEARN Machine Learning based on performance
            %
            % P = LEARN(P) The updates the aggregated agents and agent
            %   performance (B,S) using the specified learning type.
            %
            % P = LEARN(P,TYPE) will reset the learning
            %
            % References:
            % [1] Cover, T., M. (1991) Universal Portfolios
            % [2] Gyorfi, L., Udina, F., Walk, H., (2008) Experiments on universal
            %           portfolio selection using data from real markets
            % [3] Algoet, P. H., Cover, T. M., (1980) Asymptotic optimality and
            %           symptotic equipartition properties of log-optimum investments
                         
            %% input parameters
            p = varargin{1}; % (dim(x) = t) price t+1;
            % set defaults for optional inputs
            optargs = {p p.dimIntraday(3)};
            % now put these defaults into the valuesToUse cell array,
            optargs(1:nargin) = varargin;
            % Place optional args in memorable variable names
            [p, t, ret, controls, t_intra, t_day] = optargs{:};
            
            %% Machine Learning  
            % compute weights
            qH0 = p.SH_fused(t+1,:);

            %% renormalise the mixtures

            qH0 = (qH0 - mean(qH0));
            norm0 = sum(abs(qH0 - mean(qH0)));
            if norm0>eps
                qH0 = (qH0 - mean(qH0)) ./ norm0;
            else
                qH0 = zeros(size(qH0));
            end
            
            %% create the performance weighted combination of experts
            b0 = qH0*controls;
            
            %% compute normalization abs(long) + abs(short)
            nu = nansum(abs(b0));
            % renormalize controls (leverage=1)
            if nu > eps 
                b0 = (1/nu) * b0;
                % update the agent mixture weights for leverage
                qH0 = (1/nu) * qH0;
            else
               % check if active trading has commenced (past the second
               % trading day): if so then hold position from previous day
               % overnight.
%                if t <= 100
                   %update the agent mixture weights for leverage
                   qH0 = repmat(1/length(qH0),size(qH0,1),size(qH0,2));

                   %zero weights
                   b0 = zeros(size(b0));
%                else
%                   qH0 = repmat(1/length(qH0),size(qH0,1),size(qH0,2));
% 
%                   % bring previous days position over to next day
%                   b0 = p.b(t,:);
%               end
            end

           %% compute transaction costs 
           % forecast the volatiity for the first 10 timebars using
           % previous days returns using a GARCH(1,1) mode and then used
           % realised volatility method after that
           
           %get returns for last 30 time bars of previous day and compute 
           % series of returns
          if all(b0)~=0 %no transactions costs since all weights zero
            if t_intra ~= p.uniquedaytr(t_day+1) %consider the intraday and daily TC's separately
               if (t_intra - p.uniquedaytr(t_day)) < 16 %need to compute returns for day t=1 since these hav'ent been computed
        %                    % compute closing prices and returns from 30 time bars back 
        %                    % from beginning of day
        %                    Pc_intra = reshape(p.x_intraday(:,1,(p.uniquedaytr(t_day)-61):p.uniquedaytr(t_day)-1),size(p.x_intraday(:,1,(p.uniquedaytr(t_day)-61):p.uniquedaytr(t_day)-1),1),size(p.x_intraday(:,1,(p.uniquedaytr(t_day)-61):p.uniquedaytr(t_day)-1),3)); 
        %                    garchret = transpose(Pc_intra(p.liquidstocks(:,1:p.noofportfoliostocks),2:end)./Pc_intra(p.liquidstocks(:,1:p.noofportfoliostocks),1:end-1))-1;
        %                    vol = zeros(size(garchret,2),1);
        %                    for i = 1:size(garchret,2)
        %                        garchmodel = garch('GARCHLags',1,'ARCHLags',1);
        %                        [estmodel,~,~] = estimate(garchmodel, garchret(:,i),'Display','off'); %infer conditional vols from fitted garch model
        %                        [V,~] = infer(estmodel,garchret(:,i));
        %                        condvar_forecast = forecast(estmodel,1,'Y0',V);
        %                        vol(i,1) = sqrt(condvar_forecast(end)); %forecast
        %                        %sig = (log Ht - log Lt)^2/4 log 2
        %                    end
        %                else
        %                    if (t_intra - p.uniquedaytr(t_day)) < 16
                      Pc_intra = reshape(p.x_intraday(:,1,(p.uniquedaytr(t_day)-61):p.uniquedaytr(t_day)-1),size(p.x_intraday(:,1,(p.uniquedaytr(t_day)-61):p.uniquedaytr(t_day)-1),1),size(p.x_intraday(:,1,(p.uniquedaytr(t_day)-61):p.uniquedaytr(t_day)-1),3)); 
                      garchret = transpose(Pc_intra(p.liquidstocks(:,1:p.noofportfoliostocks),2:end)./Pc_intra(p.liquidstocks(:,1:p.noofportfoliostocks),1:end-1))-1;
                      vol = zeros(size(garchret,2),1);
                      for i = 1:size(garchret,2)
                        options = optimoptions(@fmincon,'Algorithm','interior-point');
                        garchmodel = garch('GARCHLags',1,'ARCHLags',1);
                        [estmodel,~,~] = estimate(garchmodel, garchret(:,i),'Display','off','options',options);
                        [V,~] = infer(estmodel,garchret(:,i)); %infer conditional vols from fitted garch model
                        condvar_forecast = forecast(estmodel,1,'Y0',V);
                        vol(i,1) = sqrt(condvar_forecast(end)); %forecast
                        %sig = (log Ht - log Lt)^2/4 log 2 %high low method
                        %using a single time bar
                      end
                else
                       Pc_intra = reshape(p.x_intraday(:,1,(p.uniquedaytr(t_day)-p.uniquedaytr(2)+1):(t_intra-p.uniquedaytr(2)+1)),...
                           size(p.x_intraday(:,1,(p.uniquedaytr(t_day)-p.uniquedaytr(2)+1):(t_intra-p.uniquedaytr(2)+1)),1),...
                           size(p.x_intraday(:,1,(p.uniquedaytr(t_day)-p.uniquedaytr(2)+1):(t_intra-p.uniquedaytr(2)+1)),3)); 
                       ret_RV = transpose(log(Pc_intra(p.liquidstocks(:,1:p.noofportfoliostocks),2:end))) - transpose(log(Pc_intra(p.liquidstocks(:,1:p.noofportfoliostocks),1:end-1)));
                      % ret_RV = ret((p.uniquedaytr(t_day)-p.uniquedaytr(2)+1):(t_intra-p.uniquedaytr(2)+1),1:end-1)';
                       vol = sum(ret_RV.^2)'; %realised volatility = sum of squared returns
               end
                   %% Intraday Transaction costs:
                   
                   % ADV - use previous 5 days to compute ADV. If 5 days have not
                   % passed then use t_days
                   if t_day <= 5
                     vol_lookback = t_day;
                     volume = reshape(p.x_intraday(p.liquidstocks(:,1:p.noofportfoliostocks),5,p.uniquedaytr(t_day-vol_lookback+1):t_intra)...
                         ,length(p.uniquedaytr(t_day-vol_lookback+1):t_intra),...
                         size(p.x_intraday(p.liquidstocks(:,1:p.noofportfoliostocks),:,:),1));
                   else
                     vol_lookback = 5; 
                     volume = reshape(p.x_intraday(p.liquidstocks(:,1:p.noofportfoliostocks),5,p.uniquedaytr(t_day-vol_lookback):t_intra)...
                         ,length(p.uniquedaytr(t_day-vol_lookback):t_intra),...
                         size(p.x_intraday(p.liquidstocks(:,1:p.noofportfoliostocks),:,:),1));               
                   end
                   ADV = sum(volume)./vol_lookback; %adv for each stocks

                   % compute number of shares traded per stock
                   % Case 1: trader - 1% of ADV per trade
                   n = (0.01/88/15).*ADV;

                   %compute transaction costs: TC = spread + sqrt(n/ADV) + indirect
                   %(fees)
                   TC = 0.007/88 + 0.0012/88 + sum(vol'*(sqrt((n.*abs(b0(1:end-1))'./ADV'))));

            else %daily trading costs
                %% Daily trading costs
                
                % Compute Transaction Cost variables
                retTC = ret;

                %specify the volatility lookback period
                if size(retTC,1) <= 15
                  vol_lookback = size(retTC,1);
                else
                   vol_lookback = 15;  
                end

                %compute volatility
                vol = std(retTC((end-vol_lookback+1):end,1:end-1),0,1);

                %ADV
                volume = reshape(p.x_daily(p.liquidstocks(:,1:p.noofportfoliostocks),5,(t_day-vol_lookback+1):t_day),vol_lookback,...
                    size(p.x_daily(p.liquidstocks(:,1:p.noofportfoliostocks),:,:),1));
                ADV = sum(volume)./vol_lookback; %adv for each stocks

               % compute number of shares traded per stock
               %% Case 1: trader - 1% of ADV per trade
               n = (0.01/88/15).*ADV;

               TC = 0.007/88 + 0.0012/88 + sum(vol.*(sqrt((n.*abs(b0(1:end-1)))./ADV)));
            end
               
            else
               TC = 0;
               n = NaN;
               vol = NaN;
           end
                           
           %% compute the leverage corrected output price relative
           % ---> update dS - portfolio wealth for current period
           if all(isnan(p.qH(t,:)))
              dS = 1;
           else
              dS = b0*(ret(end,:)'-1) + 1; % price relatives
           end

           % update the properties
           p.qH(t+1,:) = qH0;
           
           % copute returns between 4:30pm and end of day, compute the
           % associated return for this period. Then add the return for
           % closing the intraday postion to the daily return
           if t_intra == p.uniquedaytr(t_day+1)
              ret_end_of_day = [p.x_daily(p.liquidstocks(1:p.noofportfoliostocks),1,t_day)./p.x_intraday(p.liquidstocks(1:p.noofportfoliostocks),1,p.uniquedaytr(t_day+1)-1);1];
              dS = dS + -p.b(t,:)*(ret_end_of_day-1); %reverse (close i.e. -ve) the last intraday postion (time t) and add to daily return dS
           end
           
           p.b(t+1,:) = b0;
           p.S(t+1,:) = p.S(t,:) * dS;
           p.PnL(t+1,:) = dS-1;
           p.S_TC(t+1,:) = p.S_TC(t,:) * dS - TC;
           p.PnL_TC(t+1,:) = dS-1-TC;
           p.TC_vec(t+1,1) = TC;
           p.TC_vec(t+1,2) = mean(n.*abs(b0(1:end-1)));
           p.TC_vec(t+1,3) = mean(vol);
            
           %% Compute relative wealths for strategies
           ind = find(cell2mat(p.strats(:,2)),1,'last'); %find index of last strategy with 2 params
           % case for 2 param strategies
           for i = 1:ind
             p.SHave_fused(t+1,i) = mean(p.SH_fused(t+1,p.indstorekl==i)); %compute mean wealth of curreect strat for all params
           end
             
           % case for 1 param stratgeies
           for i = (ind+1):size(p.strats,1)
             p.SHave_fused(t+1,i) = mean(p.SH_fused(t+1,p.indstorel==i)); %compute mean wealth of curreect strat for all params
           end
        end
        
        %%%%%%%%%%%%%%%
        %% Partition %%
        function p = partition_Intraday(varargin)
            % LEARNING/PARTITION Partition the data
            %
            % P = PARTITION(P) single partition: [111...1] for the
            %   partition type 'trivial'.
            % 
            % P = PARTITION(P,T) Relative to time T.
            
            p = varargin{1}; 
            % set defaults for optional inputs
            optargs = {p p.dimIntraday(3)};
            % now put these defaults into the valuesToUse cell array,
            optargs(1:nargin) = varargin;
            % Place optional args in memorable variable names
            [p, t] = optargs{:};
            % create the partitions
            p.p = true(1,t);
        end

        function display(p)
            % PATTERN/DISPLAY Display the Pattern object
            %
            % See Also PATTERN 
            
            disp(p);
            if ~isempty(p.x_daily)
                fprintf('\tParameters\n');
                fprintf('\t--------------------\n');
                fprintf('\tk: %s\n',num2str(p.k));
                fprintf('\tell: %s\n',num2str(p.ell));
                fprintf('\tstrategies: %d\n',size(p.strats,1));
                fprintf('\n');
                fprintf('\tData x(M,N,P)\n');
                fprintf('\t--------------------\n');
                fprintf('\tobjects (N) : %d\n',p.dimIntraday(1));
                fprintf('\tfactors (M) : %d\n',p.dimIntraday(2));
                fprintf('\tstates  (P) : %d\n',p.dimIntraday(3));
                fprintf('\tpartitions  : %d\n',size(p.p,1));
                fprintf('\n');
            end
        end
    end
end

    %% Weight construction from signlas
    function varargout = controls_Intraday(varargin)
    % strategyfn.m calls the w-th strategy function and retruns the buy/sell signals
    % for the w-th strategy at time t. The function will also return any other
    % information required by each strategy such as the value of the strategy
    % at the previous time increment. Additinally, the function computes the
    % portfolio weights for a given expert.
    %
    % controls_Intraday(xt,strategyfn,Inputs,sigprev,ind,ell1,k1)
    %
    % controls_Intraday(xt,strategyfn,Inputs,sigprev,ind,ell1)

    %% initialise the input variables
    xt = varargin{1};
    strategyfn = varargin{2};
    Inputs = varargin{3}';
    sigprev = varargin{4};
    ind = varargin{5};
    weighttype = varargin{6};
    ell1 = varargin{7};
    optargs = {xt,strategyfn,Inputs,sigprev,ind,weighttype,ell1,{}};
    % now put these defaults into the valuesToUse cell array
    optargs(1:nargin) = varargin;

    clear varargin

    % Place optional args in memorable variable names 
    [xt,strategyfn,Inputs,sigprev,ind,weighttype,ell1,k1] = optargs{:};

    %% Compute the weights of the w-th expert/strategy

    % Find what inputs the strategy function requires
    for i = 1:length(Inputs)
       if strfind(Inputs{i},'n2')~=0 
           stratargs{1} = xt;
           stratargs{2} = ell1;
       else
           stratargs{1} = xt;
           stratargs{2} = ell1;
           stratargs{3} = k1; 
       end
    end

    % Call the w-th strategy and compute the signals
    [varargout] = strategyfn(stratargs(:));
    signals = varargout{1,1};  %n-th agent's buy/sell/hold/signals
    clear varargout stratargs

    %% Take into account previous period signals
    % Get previous period buy/sell/hold signals for expert
    signew = zeros(size(signals)); 

    if all(isnan(sigprev(ind,1:end-1))) %check if this is the case for t=tmin
      signew = signals;
    else %t>tmin
      if all(sigprev(ind,1:end-1)==0) %check if all previous signals were zeros
        signew = signals;
      else
          for i = 1:length(signals) 
            if signals(i) == 0 %check to see if current signal is a zero
              signew(i) = sign(sigprev(ind,i)); %fill in zeros for current signal with previous signs of signals
            else
              signew(i) = signals(i); 
            end
          end
      end
    end

%% Tranform Signals to Weights
%Get closing prices for stocks up to time t
Pc = reshape(xt(:,1,:),size(xt(:,1,:),1),size(xt(:,1,:),3)); %(stocks*time)

%specify the volatility lookback period
if size(Pc,2) < 90
  vol_lookback = size(Pc,2);
else
  vol_lookback = 90;  
end

%either volatility loading or inverse volatility laoding 
switch weighttype
  % Method 1: Volatility loading method
  case 'volloading'
  % weights are loaded according to stocks volatility over the last 90 days
     if all(signew==0)
       w = zeros(length(signew)+1,1); %set all weights to zero including rf asset
     elseif all(signew>=0)
       sig = std(Pc(signew>0,(end-vol_lookback+1):end),0,2);
       if sig == 0
           sig = 1;
       end
       w = [0.5*signew;-0.5]; %include rf asset
       w(signew>0) = (1/sum(sig))*w(signew>0).*sig; %normalize the weights to sum to 0
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
           error('weight not permitted')
       end
     elseif all(signew<=0)
       sig = std(Pc(signew<0,(end-vol_lookback+1):end),0,2);
       if sig == 0
           sig = 1;
       end
       w = [0.5*signew;0.5]; %include rf asset
       w(signew<0) = (1/sum(sig))*w(signew<0).*sig; %normalize the weights to sum to 0
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
           error('weight not permitted')
       end
     else
       sigpos = std(Pc(signew>0,(end-vol_lookback+1):end),0,2);
       signeg = std(Pc(signew<0,(end-vol_lookback+1):end),0,2);
       if sigpos == 0
           sigpos = 1;
       end
       if signeg == 0
           signeg = 1;
       end
       w = [signew;0]; %include rf asset
       w(signew>0) = 0.5*(1/sum(abs(sigpos)))*(sigpos).*w(signew>0);
       w(signew<0) = 0.5*(1/sum(abs(signeg)))*(signeg).*w(signew<0);
       w(end) = sum(w); %weight on rf asset doen purely to ensure weights
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
           error('weight not permitted')
       end 
     end
     
    % Method 2: Inverse Volatility loading method
    case 'invvolloading' 
     if all(signew==0)
       w = zeros(length(signew)+1,1); %set all weights to zero including rf asset
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=0
           disp('nay')
       end
     elseif all(signew>=0)
       sig = std(Pc(signew>0,(end-vol_lookback+1):end),0,2);
       w = [0.5*signew;-0.5]; %include rf asset
       w(signew>0) = (1/sum(1./sig))*w(signew>0).*(1./sig); %normalize the weights to sum to 0
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
         error('weight not permitted')
       end
     elseif all(signew<=0)
       sig = std(Pc(signew<0,(end-vol_lookback+1):end),0,2);
       w = [0.5*signew;0.5]; %include rf asset
       w(signew<0) = (1/sum(1./sig))*w(signew<0).*(1./sig); %normalize the weights to sum to 0
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
           error('weight not permitted')
       end
     else
       sigpos = std(Pc(signew>0,(end-vol_lookback+1):end),0,2);
       signeg = std(Pc(signew<0,(end-vol_lookback+1):end),0,2);
       w = [signew;0]; %include rf asset
       w(signew>0) = 0.5*(1/sum(abs(1./sigpos)))*(1./sigpos).*w(signew>0);
       w(signew<0) = 0.5*(1/sum(abs(1./signeg)))*(1./signeg).*w(signew<0);
       w(end) = sum(w); %weight on rf asset doen purely to ensure weights sum to 0
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
          error('weight not permitted')
       end
     end
end 

    varargout{1} = w;
    end   
%                           +-----+
%---------------------------| EOF |------------------------------------+
%                           +-----+