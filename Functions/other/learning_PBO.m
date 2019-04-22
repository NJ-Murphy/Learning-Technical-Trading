classdef learning_PBO
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
    % See Also: LEARNING/CONTROLS, LEARNING/OFFLINE, LEARNING/LEARN, LEARNING/ONLINE
    %
    %
    % A. OHLC patterns
    % B. Fundamental model patterns
    %
    % 1. *Data* 
    %   2.1 M  entities (objects - stocks)
    %   2.2 P  features (OHLCV)
    %   2.3 T  observations (date-times)
    % 2. *Pattern* (k-tuple) [historic or user provided]
    %   2.1 k-tuple (k free parameter)
    %       2.2.1. k=1...K
    %       2.2.2. k=1,...,N*K multiples of N (DSI)
    %   2.2 k-tuple ConSet=[A,b] A*x>=b active/absolute (0,1)
    %       2.2.1. active (many stocks,stock+cash)
    %       2.2.2. absolute (single stock, long only portfolio)
    % 4. *Distance*
    %   4.1. surface norm
    %   4.2. vector norm
    % 6. *Agents*
    %   6.1. Active
    %   6.2. Canonical Agents:
    %       6.2.1. controls H (NxMxT) M objects, N agents, T times)
    %       6.2.2. performance SH (arithmetic) (TxN for T times and N agents)
    % 7. *Algorithm*
    %   7.1. Online/Offline
    %   7.2. Parallel computing
    %   7.3. Interface agents (h,SH)
    % 8. *Machine Learning* (UNIV + active)
    %   8.1. Weighted Arithmetic Average over all agents to create optimal predictors
    %       8.1.1. Performance weighted averaging (using arithmetic averaging)
    %           B(T) = SUM(K,L) (SH(T-1|K,L) H(T|K,L)) / SUM(K,L) SH(T-1|K,L)
    %           this is probability wieghted where the probability is
    %           proportional to the returns.
    %   8.2. parameters: window W_L, forgetting factor Lamba_L
    %   8.3. either fully invested or active.
    % 9. *Sectors and States* (clusters)
    %   9.1. CxP for C clusters and M objects
    %
    %% public properties
    properties
        b = []; % aggregated controls    (T,M) Time x Objects  ---- best combined strategies 
        S = []; % aggregated performance (T,1) Time x 1  -----  weights 
        h = []; % agent controls      (N,M,T)  Agents x Objects x Time  ---- 
        SH = []; % agents performance (T,N)    Time x Agents --- Trading Strategy- time series/ strategies/ rules   
        SHave = []; %mean wealth for each strategy 
        PnL = []; %profit and losses overall strategy
        PnL_experts = []; %profit and losses experts
        S_TC = [];
        PnL_TC =[];
        TC_vec = [];
        %stockstraded = [];  %matrix to keep track of 'noofportfolio' most liquid stocks every dayinterval periods
    end
    
    %% private properties
    properties(Access = private)
        x = []; % data as price relatives (M,F,T) Objects x Features x Time
        k = []; % vector of values for parameter 2- n2
        ell = []; % vector of values for parameter 1- n/n1
        p = []; % current partition
        mnp = []; % dimensionality (M,F,T)
        weighttype = ''; %signal to weight transformation type (volatility or inverse volatility 
        qH = []; % agents controls (T,N) Time x Agents
        strats = {}; % names of the strategies
        STEFI = []; %STEFI index 
        sectors = {}; %sector clusters
        tmin0 = []; %number of periods before online learning begins
        indstorel = [];
        indstorekl = [];
        dayinterval = []; %number of days after liquidity of stocks is checked
        stocknames = []; %names of stocks needed to find clusters
        noofportfoliostocks = []; %number of stocks to be used in the portfolio construction throughout time
        liquidstocks = []; %indices of 'noofportfoliostocks' most liquid stocks at each dayinterval
        hlt = [];  %contols for strategies with 1 parameter
        hlkt = [];  %contols for strategies with 2 parameters
        noofexperts = [];
        ind_jump = [];
    end
    
    %% methods
    methods
        %% constructor
        function p = learning_PBO(varargin)
            % P = LEARNING/LEARNING constructor
            %
            % P = LEARNING(X,K,ELL,CI,NTYPE,LNTYPE,LTYPE) For data X (price relatives)
            %   and is MxNxP dimensional data for M objects, N features and
            %   P date-times. dim(X) is dim(Price)-1. A typical object is
            %   a stocks, e.g. AGL, a typical feature is a factors such as
            %   OPEN, HIGH, LOW, CLOSE (OHLC), and date-times are the time
            %   stamps for the data.
            %
            %  P = LEARNING(X,XNK,ELL,CI,NTYPE,LNTYPE,LTYPE) For matching
            %   pattern XNK instead of K-tuple range.
            %
            % The data X is homogenized in time. K is the set of tuple
            % sizes if there is no matching pattern e.g. [1:3], if there
            % is a matching pattern it is the number of matching patterns. 
            % The number of matching times is ELL this is typically 10 
            % per partition. The cluster definitions is CI, this is by 
            % default the trivial cluster (all the objects) when CI is 
            % empty. CI is KxM for K clusters of M objects. NTYPE is the 
            % strategy normalisation is 'active'. LNTYPE is the agent normalisations is 'active'.
            % The LTYPE is the learning type this is by default 'univ'.
            %
            % Note: X are price relatives (P(t)/P(t-1)). These can be
            % conveniently computed using EXP(DIFF(LOG(P))).
            %
            % See Also LEARNING/DISPLAY, LEARNING/PARTITION
            
            if nargin==0
                % null constructor
            else
                %% Compute price relatives as returns
                x0 = varargin{1}; %this corresponds to the input Data
                
                % size of the data (m=stocks, n=features, p=dates/times)
                mnp0 = size(x0);
                
                params2 = varargin{2}; %extract the parameter values for k
                params1 = varargin{3}; %extract the parameter values for ell
                
                %---> specify number of strategies to be used
                strats = varargin{4}; 
                
                %---> input STEFI
                STEFI = varargin{5};
                
                dayinterval = varargin{6};
                
                % ---> stock names
                stocknames = varargin{7};
                
                %---> sector clusters
                sectors = varargin{8};
                
                %---> stocks for portfolio contruction
                noofportfoliostocks = varargin{9};
                
                liquidstocks = varargin{10};
                
                weighttype = varargin{11};
                
                ind_jump = varargin{12};
                
                %---> set defaults for optional inputs
                optargs = {x0 params2 params1 strats STEFI dayinterval stocknames...
                    sectors noofportfoliostocks liquidstocks weighttype ind_jump};

                % Place optional args in memorable variable names
                [p.x, p.k, p.ell, p.strats, p.STEFI, p.dayinterval,...
                    p.stocknames, p.sectors, p.noofportfoliostocks,...
                    p.liquidstocks, p.weighttype, p.ind_jump] = optargs{:};

                % update the size property
                p.mnp = mnp0;
                
                %---> Number of sectors
                C = length(p.sectors);
                
                %---> control parameters
                L = size(p.ell,2); % number of param 1
                
                % number of experts 1 paramter and 2 parameters
                noofexperts1param = C*L*sum(cell2mat(strats(:,2))==0);
                noofexperts2params = C*sum(cell2mat(p.strats(:,2)))*sum(sum(p.k>max(p.ell)):sum(p.k>min(p.ell)));
                p.noofexperts = noofexperts2params + noofexperts1param;
                
                % initialise state variables and pre-allocate memory
                p.b = nan(p.mnp(3)+1,noofportfoliostocks+1); % Time x Objects (assets + rf asset)   
                p.h = nan(p.noofexperts,noofportfoliostocks+1,p.mnp(3)+1); % Agents x Objects(stocks + rf asset) x Time    
                p.qH = nan(p.mnp(3)+1,p.noofexperts); % Time x Experts
                p.S = ones(p.mnp(3)+1,1); % Time x 1
                p.S_TC = ones(p.mnp(3)+1,1); % Time x 1
                p.SH = ones(p.mnp(3)+1,p.noofexperts); %  Time x Experts  (length(p.sectors) for the stock clusters
                p.SHave = ones(p.mnp(3)+1,size(strats,1)); %Time x strategies
                %p.stockstraded = cell(noofportfoliostocks,floor(((p.mnp(3)-(max(p.k) + max(p.ell) + 1))/p.dayinterval)));
                p.hlkt = zeros(noofexperts2params,p.noofportfoliostocks+1,p.mnp(3)+1); %matrix for strats with 2 params
                p.hlt = zeros(noofexperts1param,p.noofportfoliostocks+1,p.mnp(3)+1); %matrix with strats for 1 param

                p.PnL_TC = zeros(p.mnp(3)+1,1); % Time x 1
                p.PnL = zeros(p.mnp(3)+1,1); % Time x 1
                p.TC_vec = zeros(p.mnp(3)+1,1); % Time x 1
                p.PnL_experts = zeros(p.mnp(3)+1,p.noofexperts); %  Time x Experts
              
                %% partition the data
                p = partition(p);
            end  
        end
        
        %% offline (to force offline estimation)
        function p = offline(p)
            % PATTERN/OFFLINE Offline estimation
            %
            % P = OFFLINE(P) to estimate (H,SH) for T for parameters k and
            % ell
            %
            % See Also: PATTERN/ONLINE
            
            %---> Specify the number of periods for the offline segement of
            % learning
            p.tmin0 = 2;
            T = p.mnp(3); %total number of time periods in the data
            
            if T < p.tmin0
                error('pattern:offline','Not enough Data L*K>T');
            end
            
            % Extract closing prices to compute returns
            Pc = reshape(p.x(:,1,1:p.tmin0-1),size(p.x(:,1,1:p.tmin0-1),1),size(p.x(:,1,1:p.tmin0-1),3)); %(stocks x time)

            % Compute price relatives (returns)
            ret_all = [ones(1,size(Pc,1)+1);[transpose(Pc(:,2:end)./Pc(:,1:end-1)),p.STEFI(2:p.tmin0-1,1)./p.STEFI(1:p.tmin0-1-1,1)]]; %only copute returns up to befre tmin0 as tmin0

            % implement the sequential learning for next time t
            for t = p.tmin0:T   %time loop
              Pc = reshape(p.x(:,1,t-1:t),size(p.x(:,1,t-1:t),1),size(p.x(:,1,t-1:t),3)); %closing prices time t-1 to t               
              
              % deal with the daily return jump for missing data piece
               if p.ind_jump > 1 && t == p.ind_jump
                 ret_all = [ret_all;ones(size(Pc,1)+1,1)'];
               else
                 ret_all = [ret_all; [transpose(Pc(:,end)./Pc(:,end-1)),p.STEFI(t,1)./p.STEFI(t-1,1)]];
               end
              
              p = online(p,t,ret_all);
            end    
               

        end
        
        %% online (to force online estimation)
        function p = online(varargin)
            % PATTERN/ONLINE Offline estimation
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
            [p, t, ret_all] = optargs{:};
            
            %% control parameters
            L = size(p.ell,2); % number of param 1
            W = size(p.strats,1); % number of strategies
            K = size(p.k,2); % number of params 2
            
            %% Model Identification Loop
            
            % Extract data from time t=1 to current time
            xt = p.x(:,:,1:t);
            
            % only get returns for the liquid stocks we are considering
            ret = ret_all(:,[p.liquidstocks(:,1:p.noofportfoliostocks) end]); %(time x stocks)                  
            
            if all(p.ell >= t)
                p.h(:,:,t+1) = zeros; %fill nans in with zeros as no experts trade in this case
                dSH = ones(size(p.SH(t,:)));

                %update intraday-daily wealth matrix
                p.SH(t+1,:) = p.SH(t,:).*dSH;
                p = learn(p,t,ret);
                return
            end
                    
            % Initialise the signal matrices for the t-th stime step
            C = length(p.sectors); %number of clusters
            
            %initialise indices for experts with 1 param
            hlktind = 1:length(p.sectors)*sum(cell2mat(p.strats(:,2)))*(K*(min(p.k)-min(p.ell))+(K-1)*((K)/2)-((max(p.k)-max(p.ell))-1)*(max(p.k)-max(p.ell))/2)+400;%+100; %indices for experts with 2 params
            hlktind1 = 0;  %initialise
            hltind = 0;
            expertparam2 = zeros(p.noofexperts,4);
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Expert Generation: %%
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
                            % check if enough data
                            k1 = p.k(k0); % ----> p.k will be a vector of chosen lookback periods to try for n2 for each strategy
                            if k1 > t
                                break
                            end
                            
                            %----> Condition to break out of loop if only one param
                            % is required for current strategy
                            if cell2mat(p.strats(w0,2)) == 0
                                break
                            end                       
                            k1 = p.k(k0); % ----> p.k will be a vector of chosen lookback periods to try for n2 for each strategy

                            %----> Continue loop if the k parameter is < or =
                            %to ell since ell==n/n1 cannot be > or = to k
                            if k1 > ell1
                               %Expert index KLrow(w0,k0,ell0)
                               hlktind1 = hlktind1+1;
                               Expind = hlktind(hlktind1); 
                               p.indstorekl(Expind) = w0;
                               %Compute the signals of the w-th expert for the
                               %current combination of k and ell
                               sig = controls(xt(logical(p.sectors{c}),:,:),strategyfn,p.strats{w0,3},p.hlkt(:,:,t),Expind,p.weighttype,ell1,k1); %(logical(p.sectors{c}),:,:)
                               p.hlkt(Expind,p.sectors{c}>0,t+1) = sig(1:end-1);  %Experts x stocks in cluster c
                               p.hlkt(Expind,sum(p.sectors{c}>0)+1,t+1) = sig(end); %risk free weight
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
                              sig = OnlineZAnticor(ret(:,logical([p.sectors{c},1])),size(ret,1), ell1, p.hlt(Expind,logical([p.sectors{c},1]),t));
                              p.hlt(Expind,logical([p.sectors{c}>0,1]),t+1) = sig;  %Experts x stocks in cluster c
                           elseif func2str(strategyfn) == "OnlineZBCRP"
                              sig = OnlineZBCRP(ret(:,logical([p.sectors{c},1])), ell1);
                              p.hlt(Expind,logical([p.sectors{c}>0,1]),t+1) = sig;  %Experts x stocks in cluster c
                           elseif func2str(strategyfn) == "OnlineAntiZBCRP"
                              sig = OnlineAntiZBCRP(ret(:,logical([p.sectors{c},1])), ell1);
                              p.hlt(Expind,logical([p.sectors{c}>0,1]),t+1) = sig;  %Experts x stocks in cluster c
                           else 
                              sig = controls(xt(logical(p.sectors{c}),:,:),strategyfn,p.strats{w0,3},p.hlt(:,:,t),Expind,p.weighttype,ell1);  %(logical(p.sectors{c}),:,:)
                              %aggregate current experts controsl into H
                              p.hlt(Expind,p.sectors{c}>0,t+1) = sig(1:end-1);  %Experts x stocks in cluster c
                              p.hlt(Expind,sum(p.sectors{c})+1,t+1) = sig(end); %risk free weight
                           end
                        end
                    end % ell                
                end % w          
            end % c
            
            % Compute the update performance for the prior agent step
            % ------->>> This is step 2 of the online algo (pg. 4 of
            %Gebbie et al. (2016))
            if any(isnan(p.h(:,:,t))) 
               dSH = ones(size(p.SH(t,:)));  %first iteration
            else
               dSH = (p.h(:,:,t)*(ret(end,:)'-1) + 1)';  %ret = price relatives
            end
            
            %check for case where no experts using 2 params traded during
            %current period
            if exist('expertparam2','var')== 0
                hlkt1 = zeros(size(p.h(:,:,t+1),1)-size(p.hlt(:,:,t+1),1),size(p.hlt(:,:,t+1),2));
            else
               %expertparam2 = expertparam2(1:hlktind1,:);
               hlkt1 = p.hlkt(1:hlktind1,:,t+1); %only pick out rows for which k loop was activated-this is included due            
               %to the fact that we used increments of 3 in k and ell in
               %script file
               if size([hlkt1;p.hlt(:,:,t+1)]) ~= (size(p.h(:,:,t+1),1) - size(p.hlt(:,:,t+1),1))
                   hlkt1 = [hlkt1;zeros(size(p.h(:,:,t+1),1)-size(p.hlt(:,:,t+1),1)-size(hlkt1,1),size(p.hlt(:,:,t+1),2))];
               end
            end
            
            % Update the agents at time t+1 and calculate their corresponding wealth
            p.h(:,:,t+1) = [hlkt1;p.hlt(:,:,t+1)]; %concatenate experts with 1 and 2 params resp.
            
            % Update the agent accumulated performance (geometric returns)
            %------->>> Here SH is the agent wealth- also need to take into
            % account borrowing costs so subtract STEFI when shorting risk free asset 
            p.SH(t+1,:) =  p.SH(t,:).*dSH;
            
            % profits and losses of experts
            p.PnL_experts(t+1,:) = dSH-1;

            %% online update the learning
            p = learn(p,t,ret);
        end
        
        %% universal learning
        function p = learn(varargin)
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
            %                         
            %% input parameters
            p = varargin{1}; % (dim(x) = t) price t+1;
            % set defaults for optional inputs
            optargs = {p p.mnp(3)};
            % now put these defaults into the valuesToUse cell array,
            optargs(1:nargin) = varargin;
            % Place optional args in memorable variable names
            [p, t, ret] = optargs{:};
            
            %% Machine Learning  
            % compute weights
            qH0 = p.SH(t+1,:);

            %% renormalise the mixtures

            qH0 = (qH0 - mean(qH0));
            norm0 = sum(abs(qH0 - mean(qH0)));
            if norm0>eps
                qH0 = (qH0 - mean(qH0))./norm0;
            else
                qH0 = zeros(size(qH0));
            end

            %% create the performance weighted combination of experts.
            % ONLINE
            % -----------------------------------------------------+
            b0 = qH0*p.h(:,:,t+1);
            
            %% compute normalization abs(long) + abs(short)
            nu = nansum(abs(b0));
            % renormalize controls (leverage = 1)
            if nu > eps 
                b0 = (1/nu) * b0;
                % update the agent mixture weights for leverage
                qH0 = (1/nu) * qH0;
            else
                % update the agent mixture weights for leverage
                qH0 = repmat(1/length(qH0),size(qH0,1),size(qH0,2));

                % zero weights
                b0 = zeros(size(b0));
            end
            
           if all(b0)~=0 %no transactions costs since all weights zero

            % Compute Transaction Cost variables
            retTC = ret;

            %specify the volatility lookback period
            if size(retTC,1) < 90
              vol_lookback = size(retTC,1);
            else
               vol_lookback = 90;  
            end
            
            %compute volatility
            volatility = std(retTC((end-vol_lookback+1):end,1:end-1),0,1);

            %ADV
            volume = reshape(p.x(p.liquidstocks(:,1:p.noofportfoliostocks),5,(t-vol_lookback+1):t),vol_lookback,size(p.x(p.liquidstocks(:,1:p.noofportfoliostocks),:,:),1));
            ADV = sum(volume)./vol_lookback; %adv for each stocks
            
            % No. shares traded per stock
            n = 0.00005*ADV;

            % Transaction costs
            TC = 0.0004 + 0.0001 + volatility*(sqrt(n.*abs(b0(1:end-1))'/ADV));
                
           else
               TC = 0;
           end
           
            % ---> update dS: portfolio wealth for current period (added
            % to tomorrows wealth)
            if all(isnan(p.qH(t,:)))
                dS = 1;
            else
                dS = b0*(ret(end,:)'-1)+1; % price relatives
            end
            
            % update the properties
            p.qH(t+1,:) = qH0;
            p.b(t+1,:) = b0;
            p.S_TC(t+1,:) = p.S(t,:) * (dS - TC);
            p.PnL_TC(t+1,:) = dS-1 - TC;
            p.S(t+1,:) = p.S(t,:)*dS;
            p.PnL(t+1,:) = dS-1;
            p.TC_vec(t+1,:) = TC;
            
            %% Compute relative wealths for strategies
            ind = find(cell2mat(p.strats(:,2)),1,'last'); %find index of last strategy with 2 params
            % case for 2 param strategies
            for i = 1:ind
               p.SHave(t+1,i) = mean(p.SH(t+1,find(p.indstorekl==i))); %compute mean wealth of curreect strat for all params
            end
             
            % case for 1 param stratgeies
            for i = (ind+1):size(p.strats,1)
               p.SHave(t+1,i) = mean(p.SH(t+1,find(p.indstorel==i))); %compute mean wealth of curreect strat for all params
            end
        end
        
        %%%%%%%%%%%%%%%
        %% Partition %%
        function p = partition(varargin)
            % PATTERN/PARTITION Partition the data
            %
            % P = PARTITION(P) single partition: [111...1] for the
            %   partition type 'trivial'. NP=ELL as the number of partitions.
            %   The type of partition as TYPE.
            %
            % P = PARTITION(P,T) Relative to time T.
            %
            % Note 1: Side Information based learning will partition the data based
            % on the state of the side information. The required state will
            % then be used to determine which partition to use at time T
            % for the model estimation and learning.
            
            p = varargin{1}; 
            % set defaults for optional inputs
            optargs = {p p.mnp(3)};
            % now put these defaults into the valuesToUse cell array,
            optargs(1:nargin) = varargin;
            % Place optional args in memorable variable names
            [p, t] = optargs{:};
            
            % create the partitions
            p.p = true(1,t);
        end

        function display(p)           
            disp(p);
            if ~isempty(p.x)
                fprintf('\tParameters\n');
                fprintf('\t--------------------\n');
                fprintf('\tk-tuples       : %s\n',num2str(p.k));
                fprintf('\tell neighbours : %s\n',num2str(p.ell));
                fprintf('\tstrategies      : %d\n',size(p.strats,1));
                fprintf('\n');
                fprintf('\tData x(M,N,P)\n');
                fprintf('\t--------------------\n');
                fprintf('\tobjects (N) : %d\n',p.mnp(1));
                fprintf('\tfactors (M) : %d\n',p.mnp(2));
                fprintf('\tstates  (P) : %d\n',p.mnp(3));
                fprintf('\tpartitions  : %d\n',size(p.p,1));
                fprintf('\n');
            end
        end
    end
end

    %% Weight construction from signlas
    function varargout = controls(varargin)
    % controls.m calls the w-th strategy function and retruns the buy/sell signals
    % for the w-th strategy at time t. The function will also return any other
    % information required by each strategy such as the value of the strategy
    % at the previous time increment. Additinally, the function computes the
    % portfolio weights for a given expert.
    %
    % controls(xt,strategyfn,Inputs,sigprev,ind,ell1,k1)
    %
    % controls(xt,strategyfn,Inputs,sigprev,ind,ell1)

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
if size(Pc,2) < 120
  vol_lookback = size(Pc,2);
else
  vol_lookback = 120;  
end

%either volatility loading or inverse volatility laoding 
switch weighttype
  % Method 2: Volatility loading method
  case 'volloading'
  % weights are loaded according to stocks volatility over the last 90 days
     if all(signew==0)
       w = zeros(length(signew)+1,1); %set all weights to zero including rf asset
     elseif all(signew>=0)
       sig = std(Pc(signew>0,(end-vol_lookback+1):end),0,2);
       w = [0.5*signew;-0.5]; %include rf asset
       w(signew>0) = (1/sum(sig))*w(signew>0).*sig; %normalize the weights to sum to 0
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
           error('weight not permitted')
       end
     elseif all(signew<=0)
       sig = std(Pc(signew<0,(end-vol_lookback+1):end),0,2);
       w = [0.5*signew;0.5]; %include rf asset
       w(signew<0) = (1/sum(sig))*w(signew<0).*sig; %normalize the weights to sum to 0
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
           error('weight not permitted')
       end
     else
       sigpos = std(Pc(signew>0,(end-vol_lookback+1):end),0,2);
       signeg = std(Pc(signew<0,(end-vol_lookback+1):end),0,2);
       w = [signew;0]; %include rf asset
       w(signew>0) = 0.5*(1/sum(abs(sigpos)))*(sigpos).*w(signew>0);
       w(signew<0) = 0.5*(1/sum(abs(signeg)))*(signeg).*w(signew<0);
       w(end) = sum(w); %weight on rf asset doen purely to ensure weights
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=1
           error('weight not permitted')
       end 
     end
     
    % Method 3: Inverse Volatility loading method
    case 'invvolloading' 
     if all(signew==0)
       w = zeros(length(signew)+1,1); %set all weights to zero including rf asset
       if round(sum(w),4)~= 0 || round(sum(abs(w)),4)~=0
           disp('weight not permitted')
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
       %sig = std(Pc(signew~=0,(end-vol_lookback+1):end),0,2);
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
%--------------------------------------------------------------------+