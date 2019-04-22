function vout = indicators(vin,mode,varargin)
%INDICATORS calculates various technical indicators
% 
% Description
%     INDICATORS is a technical analysis tool that calculates various
%     technical indicators.  Technical analysis is the forecasting of
%     future financial price movements based on an examination of past
%     price movements.  Most technical indicators require at least 1
%     variable argument.  If these arguments are not supplied, default
%     values are used.
% 
% Syntax
%     Momentum
%         cci                  = indicators([hi,lo,cl]      ,'cci'    ,tp_per,md_per,const)
%         roc                  = indicators(price           ,'roc'    ,period)
%         rsi                  = indicators(price           ,'rsi'    ,period)
%         [fpctk,fpctd]        = indicators([hi,lo,cl]      ,'fsto'   ,k,d)
%         [spctk,spctd]        = indicators([hi,lo,cl]      ,'ssto'   ,k,d)
%         [fpctk,fpctd,jline]  = indicators([hi,lo,cl]      ,'kdj'    ,k,d)
%         willr                = indicators([hi,lo,cl]      ,'william',period)
%         [dn,up,os]           = indicators([hi,lo]         ,'aroon'  ,period)
%         tsi                  = indicators(cl              ,'tsi'    ,r,s)
%     Trend
%         sma                  = indicators(price           ,'sma'    ,period)
%         ema                  = indicators(price           ,'ema'    ,period)
%         [macd,signal,macdh]  = indicators(cl              ,'macd'   ,short,long,signal)
%         [pdi,mdi,adx]        = indicators([hi,lo,cl]      ,'adx'    ,period)
%         t3                   = indicators(price           ,'t3'     ,period,volfact)
%     Volume
%         obv                  = indicators([cl,vo]         ,'obv')
%         cmf                  = indicators([hi,lo,cl,vo]   ,'cmf'    ,period)
%         force                = indicators([cl,vo]         ,'force'  ,period)
%         mfi                  = indicators([hi,lo,cl,vo]   ,'mfi'    ,period)
%     Volatility
%         [middle,upper,lower] = indicators(price           ,'boll'   ,period,weight,nstd)
%         [middle,upper,lower] = indicators([hi,lo,cl]      ,'keltner',emaper,atrmul,atrper)
%         atr                  = indicators([hi,lo,cl]      ,'atr'    ,period)
%         vr                   = indicators([hi,lo,cl]      ,'vr'     ,period)
%         hhll                 = indicators([hi,lo]         ,'hhll'   ,period)
%     Other
%         [index,value]        = indicators(price           ,'zigzag' ,moveper)
%         change               = indicators(price           ,'compare')
%         [pivot sprt res]     = indicators([dt,op,hi,lo,cl],'pivot'  ,type)
%         sar                  = indicators([hi,lo]         ,'sar'    ,step,maximum)
% 
% Arguments
%     Outputs
%         cci/roc/rsi/willr/tsi/sma/ema/t3/obv/cmf/force/mfi/atr/change/sar/vr/hhll
%                 - single output vector
%         macd/signal/macdh
%                 - moving average convergence divergence vector/
%                   signal line/ macd histogram
%         middle/upper/lower
%                 - middle/upper/lower band for bollinger bands or keltner
%                   channels
%         fpctk/fpctd/spctk/spctd/jline
%                 - fast/slow percent k/d and the J Line
%         index/value
%                 - index and value for each point in zigzag
%         pivot/sprt/res
%                 - vectors for pivot point, support lines, resistance
%                 lines
%         dn/up/os
%                 - aroon down/aroon up/aroon oscillator
%     Inputs
%         price   - any price vector (e.g. open, high, ...)
%         dt/op/hi/lo/cl/vo
%                 - matlab serial date, open/high/low/close price, and volume of data
%     Mode
%         Momentum
%             cci     - Commodity Channel Index
%             roc     - Rate of Change
%             rsi     - Relative Strength Index
%             fsto    - Fast Stochastic Oscillator
%             ssto    - Slow Stochastic Oscillator
%             kdj     - KDJ Indicator
%             william - William's %R
%             aroon   - Aroon
%             tsi     - True Strength Index
%         Trend
%             sma     - Simple Moving Average
%             ema     - Exponential Moving Average
%             macd    - Moving Average Convergence Divergence
%             adx     - Wildmer's DMI (ADX)
%             t3      - Triple EMA (Not the same as EMA3)
%         Volume
%             obv     - On-Balance Volume
%             cmf     - Chaikin Money Flow
%             force   - Force Index
%             mfi     - Money Flow Index
%         Volatility
%             boll    - Bollinger Bands
%             keltner - Keltner Channels
%             atr     - Average True Range
%             vr      - Volatility Ratio
%             hhll    - Highest High, Lowest Low
%         Other
%             zigzag  - ZigZag
%             compare - relative price compared to first input
%             pivot   - Pivot Points
%             sar     - Parabolic SAR (Stop And Reverse)
% 
%     Variable Arguments
%         period  - number of periods over which to make calculations
%         short/long/signal
%                 - number of periods for short/long/signal for macd
%         period/weight/nstd
%                 - number of periods/weight factor/number of standard
%                   deviations
%         k/d/j   - number of periods for %K/%D/JLine
%         r/s     - number of periods for momentum/smoothed momentum
%         emaper/atrmul/atrper
%                 - number of periods for ema/atr multiplier/number of
%                   periods for atr for keltner
%         tp_per/md_per/const
%                 - number of periods for true price/number of periods for
%                   mean deviation/constant for cci
%         moveper - movement percent for zigzag
%         type    - string for pivot point method pick one of 's', 'f', 'd'
%                   which stand for 'standard', 'fibonacci', 'demark'
%         step/maximum
%                 - value to add to acceleration factor at each increment
%                   maximum value acceleration factor can reach
%         volfact - volume factor for t3
% 
% Notes
%     - there are no argument checks
%     - all prices must be column oriented
%     - the parabolic sar indicator is not completely correct
%     - if there is a tie between price points, the aroon indicator uses
%     the one farthest back
%     - 2 methods are available to calculate the ema for the tsi
%     indicator.  Simply uncomment whichever one is desired.
%     - the t3 indicator uses a different ema than ta-lib
%     - 3 methods are available to calculate the kdj.  Simply uncomment
%     whichever one is desired.
% 
% Example
%     load disney.mat
%     vout = indicators([dis_HIGH,dis_LOW,dis_CLOSE],'fsto',14,3);
%     fpctk = vout(:,1);
%     fpctd = vout(:,2);
%     plot(1:length(fpctk),fpctk,'b',1:length(fpctd),fpctd,'g')
%     title('Fast Stochastics for Disney')
% 
% Further Information
%     For an in depth analysis of how many of these technical indicators
%     work, refer to one of the following websites or Google it.
%     http://stockcharts.com/
%     http://www.investopedia.com/
%     http://www.ta-lib.org/
%

% Version : 1.1.3 (05/24/2013)
% Author  : Nate Jensen
% Created : 10/10/2011
% History :
%  - v1.0 10/25/2011 : initial release of 21 indicators
%  - v1.1 03/04/2012 : 23 indicators, fixed date conversion issue
%  - v1.1.1 03/25/2012 : 24 indicators
%  - v1.1.2 03/21/2013 : 25 indicators
%  - v1.1.3 05/24/2013 : 27 indicators, bug fixes

% To Do
%  - add more indicators

%%% Main Function

% Initialize output vector
vout = [];

% Number of observations
observ = size(vin,1);

% Switch between the various modes
switch lower(mode)
    
    %%% Momentum
%==========================================================================
    case 'cci'      % Commodity Channel Index
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable Argument Input
        if isempty(varargin)
            tp_per = 20;
            md_per = 20;
            const  = 0.015;
        else
            tp_per = varargin{1};
            md_per = varargin{2};
            const  = varargin{3};
        end
        
        % Typical Price
        tp = (hi+lo+cl)/3;
        
        % Simple moving average of typical price
        smatp = sma(tp,tp_per,observ);
        
        % Sum of the mean absolute deviation
        smad = nan(observ,1);
        cci  = smad;    % preallocate cci
        for i1 = md_per:observ
            smad(i1) = sum(abs(smatp(i1)-tp(i1-md_per+1:i1)));
        end
        
        % Commodity Channel Index
        i1 = md_per:observ;
        cci(i1) = (tp(i1)-smatp(i1))./(const*smad(i1)/md_per);
        
        % Format Output
        vout = cci;
        
    case 'roc'      % Rate of Change
        % Input Data
        cl = vin;
        
        % Variable Argument Input
        if isempty(varargin)
            period = 12;
        else
            period = varargin{1};
        end
        
        % Rate of Change
        roc = nan(observ,1);
        % calculate rate of change
        roc(period+1:observ) = ((cl(period+1:observ)- ...
            cl(1:observ-period))./cl(1:observ-period))*100;
        
        % Format Output
        vout = roc;
        
    case 'rsi'      % Relative Strength Index
        % Input Data
        cl = vin;
        
        % Variable Argument Input
        if isempty(varargin)
            period = 14;
        else
            period = varargin{1};
        end
        
        % Determine how many nans are in the beginning
        nanVals  = isnan(cl);
        firstVal = find(nanVals == 0, 1, 'first');
        numLeadNans = firstVal - 1;
        
        % Create vector of non-nan closing prices
        nnanvin = cl(~isnan(cl));
        
        % Take a diff of the non-nan closing prices
        diffdata    = diff(nnanvin);
        priceChange = abs(diffdata);
        
        % Create '+' Delta vectors and '-' Delta vectors
        advances = priceChange;
        declines = priceChange;
        
        advances(diffdata < 0)  = 0;
        declines(diffdata >= 0) = 0;
        
        % Calculate the RSI of the non-nan closing prices. Ignore first non-nan
        % vin b/c it is a reference point. Take into account any leading nans
        % that may exist in vin vector.
        trsi = nan(size(diffdata, 1)-numLeadNans, 1);
        for i1 = period:size(diffdata, 1)
            % Gains/losses
            totalGain = sum(advances((i1 - (period-1)):i1));
            totalLoss = sum(declines((i1 - (period-1)):i1));

            % Calculate RSI
            rs         = totalGain ./ totalLoss;
            trsi(i1) = 100 - (100 / (1+rs));
        end
        
        % Pre allocate vector taking into account reference value and leading nans.
        % length of vector = length(vin) - # of reference values - # of leading nans
        rsi = nan(size(cl, 1)-1-numLeadNans, 1);
        
        % Populate RSI
        rsi(~isnan(cl(2+numLeadNans:end))) = trsi;
        
        % Format Output
        vout = [nan(numLeadNans+1, 1); rsi];
        
    case 'fsto'     % Fast Stochastic Oscillator
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable Argument Input
        if isempty(varargin)
            kperiods = 14;
            dperiods = 3;
        else
            kperiods = varargin{1};
            dperiods = varargin{2};
        end
        
        % Fast %K
        fpctk = nan(observ,1);                      % preallocate Fast %K
        llv = zeros(observ,1);                      % preallocate lowest low
        llv(1:kperiods) = min(lo(1:kperiods));      % lowest low of first kperiods
        for i1 = kperiods:observ                    % cycle through rest of data
            llv(i1) = min(lo(i1-kperiods+1:i1));    % lowest low of previous kperiods
        end
        hhv = zeros(observ,1);                      % preallocate highest high
        hhv(1:kperiods) = max(hi(1:kperiods));      % highest high of first kperiods
        for i1 = kperiods:observ                    % cycle through rest of data
            hhv(i1) = max(hi(i1-kperiods+1:i1));    % highest high of previous kperiods
        end
        nzero        = find((hhv-llv) ~= 0);
        fpctk(nzero) = ((cl(nzero)-llv(nzero))./(hhv(nzero)-llv(nzero)))*100;
        
        % Fast %D
        fpctd                = nan(size(cl));
        fpctd(~isnan(fpctk)) = ema(fpctk(~isnan(fpctk)),dperiods, ...
            length(fpctk(~isnan(fpctk))));
        
        % Format Output
        vout = [fpctk,fpctd];
        
    case 'ssto'     % Slow Stochastic Oscillator
        % Input data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable Argument Input
        if isempty(varargin)
            kperiods = 14;
            dperiods = 3;
        else
            kperiods = varargin{1};
            dperiods = varargin{2};
        end
        
        % Fast %K
        fpctk = nan(observ,1);                      % preallocate Fast %K
        llv = zeros(observ,1);                      % preallocate lowest low
        llv(1:kperiods) = min(lo(1:kperiods));     % lowest low of first kperiods
        for i1 = kperiods:observ                    % cycle through rest of data
            llv(i1) = min(lo(i1-kperiods+1:i1));   % lowest low of previous kperiods
        end
        hhv = zeros(observ,1);                      % preallocate highest high
        hhv(1:kperiods) = max(hi(1:kperiods));    % highest high of first kperiods
        for i1 = kperiods:observ                    % cycle through rest of data
            hhv(i1) = max(hi(i1-kperiods+1:i1));  % highest high of previous kperiods
        end
        nzero        = find((hhv-llv) ~= 0);
        fpctk(nzero) = ((cl(nzero)-llv(nzero))./(hhv(nzero)-llv(nzero)))*100;
        
        % Fast %D
        fpctd                = nan(size(cl));
        fpctd(~isnan(fpctk)) = ema(fpctk(~isnan(fpctk)),dperiods, ...
            length(fpctk(~isnan(fpctk))));
        
        % Slow %K
        spctk = fpctd;
        
        % Slow %D
        spctd = nan(size(cl));
        spctd(~isnan(spctk)) = ema(spctk(~isnan(spctk)),dperiods, ...
            length(spctk(~isnan(spctk))));
        
        % Format Output
        vout = [spctk,spctd];
        
    case 'kdj'     % KDJ Indicator
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable Argument Input
        if isempty(varargin)
            kperiods = 14;
            dperiods = 3;
            % jperiods = 5;
        else
            kperiods = varargin{1};
            dperiods = varargin{2};
            % jperiods = varargin{3};
        end
        
        % Fast %K
        fpctk = nan(observ,1);                      % preallocate Fast %K
        llv = zeros(observ,1);                      % preallocate lowest low
        llv(1:kperiods) = min(lo(1:kperiods));     % lowest low of first kperiods
        for i1 = kperiods:observ                    % cycle through rest of data
            llv(i1) = min(lo(i1-kperiods+1:i1));   % lowest low of previous kperiods
        end
        hhv = zeros(observ,1);                      % preallocate highest high
        hhv(1:kperiods) = max(hi(1:kperiods));    % highest high of first kperiods
        for i1 = kperiods:observ                    % cycle through rest of data
            hhv(i1) = max(hi(i1-kperiods+1:i1));  % highest high of previous kperiods
        end
        nzero        = find((hhv-llv) ~= 0);
        fpctk(nzero) = ((cl(nzero)-llv(nzero))./(hhv(nzero)-llv(nzero)))*100;
        
        % Fast %D
        fpctd                = nan(size(cl));
        fpctd(~isnan(fpctk)) = ema(fpctk(~isnan(fpctk)),dperiods,observ);
        
        % Method # 1:
        jline = 3*fpctk-2*fpctd;
        
        % Method # 2:
        % jline                = nan(size(cl));
        % jline(~isnan(fpctk)) = ema(fpctk(~isnan(fpctk)),jperiods,observ);
        
        % Method # 3:
        % Slow %K
        % spctk = fpctd;
        
        % Slow %D
        % spctd = nan(size(cl));
        % spctd(~isnan(spctk)) = ema(spctk(~isnan(spctk)),dperiods,observ-dperiods+1);
        
        % J Line
        % jline = 3*spctk-2*spctd;
        
        % Format Output
        % vout = [spctk,spctd,jline];
        
        % Format Output
        vout = [fpctk,fpctd,jline];
        
    case 'william'  % William's %R
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable Argument Input
        if isempty(varargin)
            period = 14;
        else
            period = varargin{1};
        end
        
        % Highest High and Lowest Low
        llv = zeros(observ,1);                      % preallocate lowest low
        llv(1:period) = min(lo(1:period));     % lowest low of first kperiods
        for i1 = period:observ                    % cycle through rest of data
            llv(i1) = min(lo(i1-period+1:i1));   % lowest low of previous kperiods
        end
        hhv = zeros(observ,1);                      % preallocate highest high
        hhv(1:period) = max(hi(1:period));    % highest high of first kperiods
        for i1 = period:observ                    % cycle through rest of data
            hhv(i1) = max(hi(i1-period+1:i1));  % highest high of previous kperiods
        end
        
        % Williams %R
        wpctr        = nan(observ,1);
        nzero        = find((hhv-llv) ~= 0);
        wpctr(nzero) = ((hhv(nzero)-cl(nzero))./(hhv(nzero)-llv(nzero))) * -100;
        
        % Format output
        vout = wpctr;
        
    case 'aroon'    % Aroon
        % Input Data
        hi     = vin(:,1);
        lo      = vin(:,2);
        
        % Variable Argument Input
        if isempty(varargin)
            period = 25;
        else
            period = varargin{1};
        end
        
        % Cumulative sum of end indices
        % Output looks like:
        % [1 16 31 46 61 76 91 ... ]
        temp_var1 = cumsum([1;(period+1:observ)'-(1:observ-period)'+1]);
        % Vector of moving indices
        % Output looks like:
        % [1 2 3 4 5 2 3 4 5 6 3 4 5 6 7 4 5 6 7 8 ... ]
        temp_var2 = ones(temp_var1(observ-period+1)-1,1);
        temp_var2(temp_var1(1:observ-period)) = 1-period;
        temp_var2(1) = 1;
        temp_var2 = cumsum(temp_var2);
        
        % Days since last n periods high/low
        [~,min_idx] = min(lo(reshape(temp_var2,period+1,observ-period)),[],1);
        [~,max_idx] = max(hi(reshape(temp_var2,period+1,observ-period)),[],1);
        
        % Aroon Down/Up/Oscillator
        aroon_dn = [nan(period,1); ((period-(period+1-min_idx'))/period)*100];
        aroon_up = [nan(period,1); ((period-(period+1-max_idx'))/period)*100];
        aroon_os = aroon_up-aroon_dn;
        
        % Format Output
        vout = [aroon_dn,aroon_up,aroon_os];
        
    case 'tsi' % True Strength Index
        % Input Data
        cl = vin(:,1);
        
        % Variable Argument Input
        if isempty(varargin)
            slow = 25;
            fast = 13;
        else
            slow = varargin{1};
            fast = varargin{2};
        end
        
        % If the lag is greater than or equal to the number of observations
        if slow >= observ || fast >= observ
            return
        end
        
        % Momentum
        mtm    = [0; (cl(2:end,1)) - cl(1:end-1,1)];
        absmtm = abs(mtm);
        
        % Calculate the exponential percentage
        k1 = 2/(slow+1);
        k2 = 2/(fast+1);
        
        % Wikipedia method for calculating ema
        % Preallocate
        ema1 = zeros(observ,1);
        ema2 = ema1;
        ema3 = ema1;
        ema4 = ema1;
        
        % EMA's
        for i1 = 2:observ
            ema1(i1) = k1 * (mtm(i1)-ema1(i1-1))    + ema1(i1-1);
            ema2(i1) = k2 * (ema1(i1)-ema2(i1-1))   + ema2(i1-1);
            ema3(i1) = k1 * (absmtm(i1)-ema3(i1-1)) + ema3(i1-1);
            ema4(i1) = k2 * (ema3(i1)-ema4(i1-1))   + ema4(i1-1);
        end
        
        % True Strength Index
        tsi = 100*ema2./ema4;
        
%         % Matlab method for calculating ema
%         % Preallocate EMA's
%         ema1 = nan(observ,1);
%         ema2 = ema1;
%         ema3 = ema1;
%         ema4 = ema1;
%         
%         % Calculate the simple moving average for the first 'exp mov avg' value.
%         ema1(slow) = sum(mtm(1:slow))/slow;
%         ema3(slow) = sum(absmtm(1:slow))/slow;
%         
%         % K*vin; 1-k
%         kvin1 = mtm(slow:observ) * k1;
%         oneK1 = 1-k1;
%         kvin3 = absmtm(slow:observ) * k1;
%         oneK2 = 1-k2;
%         
%         % First period calculation
%         ema1(slow) = kvin1(1) + (ema1(slow) * oneK1);
%         ema3(slow) = kvin3(1) + (ema3(slow) * oneK1);
%         
%         % Remaining periods calculation
%         for i1 = slow+1:observ
%             ema1(i1) = kvin1(i1-slow+1) + (ema1(i1-1) * oneK1);
%             ema3(i1) = kvin3(i1-slow+1) + (ema3(i1-1) * oneK1);
%         end
%         
%         % Calculate the simple moving average for the first 'exp mov avg' value.
%         ema2(slow+fast-1) = sum(ema1(slow:slow+fast-1))/fast;
%         ema4(slow+fast-1) = sum(ema3(slow:slow+fast-1))/fast;
%         
%         % K*vin; 1-k
%         kvin2 = ema1(slow+fast-1:observ) * k2;
%         kvin4 = ema3(slow+fast-1:observ) * k2;
%         
%         % First period calculation
%         ema2(slow+fast-1) = kvin2(1) + (ema2(slow+fast-1) * oneK2);
%         ema4(slow+fast-1) = kvin4(1) + (ema4(slow+fast-1) * oneK2);
%         
%         % Remaining periods calculation
%         for i1 = slow+fast:observ
%             ema2(i1) = kvin2(i1-fast-slow+2) + (ema2(i1-1) * oneK2);
%             ema4(i1) = kvin4(i1-fast-slow+2) + (ema4(i1-1) * oneK2);
%         end
%         
%         % True Strength Index
%         tsi = 100*ema2./ema4;
        
        % Format Output
        vout = tsi;
        
%--------------------------------------------------------------------------
    
    %%% Trend
%==========================================================================
    case 'sma'      % Simple Moving Average
        % Input Data
        price = vin;
        
        % Variable Argument Input
        if isempty(varargin)
            period = 20;
        else
            period = varargin{1};
        end
        
        % Simple Moving Average
        simmovavg = sma(price,period,observ);
        
        % Format Output
        vout = simmovavg;
        
    case 'ema'      % Exponential Moving Average
        % Input Data
        price = vin;
        
        % Variable Argument Input
        if isempty(varargin)
            period = 20;
        else
            period = varargin{1};
        end
        
        % Exponential Moving Average
        expmovavg = ema(price,period,observ);
        
        % Format Output
        vout = expmovavg;
        
    case 'macd'     % Moving Average Convergence Divergence
        % Input Data
        cl  = vin;
        
        % Variable Argument Input
        if isempty(varargin)
            short  = 12;
            long   = 26;
            signal = 9;
        else
            short  = varargin{1};
            long   = varargin{2};
            signal = varargin{3};
        end
        
        % EMA of Long Period
        [ema_lp status] = ema(cl,long,observ);
        if ~status
            return
        end
        
        % EMA of Short Period
        ema_sp = ema(cl,short,observ);
        
        % MACD
        MACD = ema_sp-ema_lp;
        
        % Signal
        [signal status] = ema(MACD(~isnan(MACD)),signal,observ-long+1);
        if ~status
            return
        end
        signal = [nan(long-1,1);signal];
        
        % MACD Histogram
        MACD_h = MACD-signal;
        
        % Format Output
        vout = [MACD,signal,MACD_h];
        
    case 'adx'      % Wilder's DMI (ADX)
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable argument input
        if isempty(varargin)
            period = 14;
        else
            period = varargin{1};
        end
        
        % True range
        h_m_l = hi-lo;                                   % high - low
        h_m_c = [0;abs(hi(2:observ)-cl(1:observ-1))];  % abs(high - close)
        l_m_c = [0;abs(lo(2:observ)-cl(1:observ-1))];   % abs(low - close)
        tr = max([h_m_l,h_m_c,l_m_c],[],2);                 % true range
        
        % Directional Movement
        h_m_h = hi(2:observ)-hi(1:observ-1);            % high - high
        l_m_l = lo(1:observ-1)-lo(2:observ);              % low - low
        pdm1  = zeros(observ-1,1);                          % preallocate pdm1
        max_h = max(h_m_h,0);
        pdm1(h_m_h > l_m_l) = max_h(h_m_h > l_m_l);         % plus
        mdm1  = zeros(observ-1,1);                          % preallocate mdm1
        max_l = max(l_m_l,0);
        mdm1(l_m_l > h_m_h) = max_l(l_m_l > h_m_h);         % minus
        pdm1 = [nan;pdm1];
        mdm1 = [nan;mdm1];
        
        % Preallocate 14 period tr, pdm, mdm, adx
        tr14  = nan(observ,1);  % 14 period true range
        pdm14 = tr14;           % 14 period plus directional movement
        mdm14 = tr14;           % 14 period minus directional movement
        adx   = tr14;           % average directional index
        
        % Calculate tr14, pdm14, mdm14, pdi14, mdi14, dmx
        tr14(period+1)  = sum(tr(period+1-period+1:period+1));
        pdm14(period+1) = sum(pdm1(period+1-period+1:period+1));
        mdm14(period+1) = sum(mdm1(period+1-period+1:period+1));
        for i1 = period+2:observ
            tr14(i1)  = tr14(i1-1)-tr14(i1-1)/period+tr(i1);
            pdm14(i1) = pdm14(i1-1)-pdm14(i1-1)/period+pdm1(i1);
            mdm14(i1) = mdm14(i1-1)-mdm14(i1-1)/period+mdm1(i1);
        end
        pdi14 = 100*pdm14./tr14;                    % 14 period plus directional indicator
        mdi14 = 100*mdm14./tr14;                    % 14 period minus directional indicator
        dmx   = 100*abs(pdi14-mdi14)./(pdi14+mdi14);% directional movement index
        
        % Average Directional Index
        adx(2*period) = sum(dmx(period+1:2*period))/(2*period-period-1);
        for i1 = 2*period+1:observ
            adx(i1) = (adx(i1-1)*(period-1)+dmx(i1))/period;
        end
        
        % Format Output
        vout = [pdi14,mdi14,adx];
        
    case 't3'       % T3
        % Input Data
        price = vin;
        
        % Variable Argument Input
        if isempty(varargin)
            period  = 5;
            volfact = 0.7;
        else
            period   = varargin{1};
            volfact = varargin{2};
        end
        
        % EMA
        ema1 = ema(price,period,observ);
        ema2 = [nan(period,1); ema(ema1(~isnan(ema1)),period,observ-period+1)];
        ema3 = [nan(2*period-1,1); ema(ema2(~isnan(ema2)),period,observ-2*period+1)];
        ema4 = [nan(3*period-1,1); ema(ema3(~isnan(ema3)),period,observ-3*period+1)];
        ema5 = [nan(4*period-1,1); ema(ema4(~isnan(ema4)),period,observ-4*period+1)];
        ema6 = [nan(5*period-1,1); ema(ema5(~isnan(ema5)),period,observ-5*period+1)];
        
        % Constants
        c1 = -(volfact*volfact*volfact);
        c2 = 3*(volfact*volfact-c1);
        c3 = -6*volfact*volfact-3*(volfact-c1);
        c4 = 1+3*volfact-c1+3*volfact*volfact;
        
        % T3
        t3 = c1*ema6+c2*ema5+c3*ema4+c4*ema3;
        
        % Format Output
        vout = t3;
        
%--------------------------------------------------------------------------
        
    %%% Volume
%==========================================================================
    case 'obv'      % On-Balance Volume
        % Input data
        cl = vin(:,1);
        vo = vin(:,2);
        
        % On-Balance Volume
        obv = vo;
        for i1 = 2:observ
            if     cl(i1) > cl(i1-1)
                obv(i1) = obv(i1-1)+vo(i1);
            elseif cl(i1) < cl(i1-1)
                obv(i1) = obv(i1-1)-vo(i1);
            elseif cl(i1) == cl(i1-1)
                obv(i1) = obv(i1-1);
            end
        end
        
        % Format Output
        vout = obv;
        
    case 'cmf'      % Chaikin Money Flow
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        vo = vin(:,4);
        
        % Variable Argument Input
        if isempty(varargin)
            period = 20;
        else
            period = varargin{1};
        end
        
        % Money Flow Multiplier
        mfm = ((cl-lo)-(hi-cl))/(hi-lo);
        
        % Money Flow Volume
        mfv = mfm*vo;
        
        % Chaikin Money Flow
        cmf = nan(observ,1);
        for i1 = period:observ
            cmf(i1) = sum(mfv(i1-period+1:i1))/sum(vo(i1-period+1:i1));
        end
        
        % Format Output
        vout = cmf;
        
    case 'force'    % Force Index
        % Input Data
        cl = vin(:,1);
        vo = vin(:,2);
        
        % Variable Argument Input
        if isempty(varargin)
            period = 13;
        else
            period = varargin{1};
        end
        
        % Force Index
        force = [nan; (cl(2:observ)-cl(1:observ-1)).*vo(2:observ)];
        force = [nan; ema(force(2:observ),period,observ-1)];
        
        % Format Output
        vout = force;
        
    case 'mfi'      % Money Flow Index
        % Input Data
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        vo = vin(:,4);
        
        if isempty(varargin)
            period = 14;
        else
            period = varargin{1};
        end
        
        % Typical Price
        tp = (hi+lo+cl)/3;
        
        % Up or Down
        upordn = ones(observ-1,1);
        upordn(tp(2:observ) <= tp(1:observ-1)) = -1;
        
        % Raw Money Flow
        rmf = tp(2:observ).*vo(2:observ);
        
        % Positive Money Flow
        pmf = zeros(observ-1,1);
        pmf(upordn == 1) = rmf(upordn == 1);
        
        % Negative Money Flow
        nmf = zeros(observ-1,1);
        nmf(upordn == -1) = rmf(upordn == -1);
        
        % Cumulative sum of end indices
        % Output looks like:
        % [1 16 31 46 61 76 91 ... ]
        temp_var1 = cumsum([1;(period:observ-1)'-(1:observ-period)'+1]);
        % Vector of moving indices
        % Output looks like:
        % [1 2 3 4 5 2 3 4 5 6 3 4 5 6 7 4 5 6 7 8 ... ]
        temp_var2 = ones(temp_var1(observ-period+1)-1,1);
        temp_var2(temp_var1(1:observ-period)) = 2-period;
        temp_var2(1) = 1;
        temp_var2 = cumsum(temp_var2);
        
        % Money Flow Ratio
        mfr = sum(pmf(reshape(temp_var2,period,observ-period)),1)'./ ...
            sum(nmf(reshape(temp_var2,period,observ-period)),1)';
        mfr = [nan(period,1); mfr];
        
        % Money Flow Index
        mfi = 100-100./(1+mfr);
        
        % Format Output
        vout = mfi;
        
%--------------------------------------------------------------------------
        
    %%% Volatility
%==========================================================================
    case 'boll'     % Bollinger Bands
        % Input data
        cl = vin;
        
        % Variable argument input
        if isempty(varargin)
            period = 20;
            weight = 0;
            nstd   = 2;
        else
            period = varargin{1};
            weight = varargin{2};
            nstd   = varargin{3};
        end
        
        % Create output vectors.
        mid  = nan(size(cl, 1), 1);
        uppr = mid;
        lowr = mid;
        
        % Create weight vector.
        wtsvec = ((1:period).^weight) ./ (sum((1:period).^weight));
        
        % Save the original data and remove NaN's from the data to be processed.
        nnandata = cl(~isnan(cl));
        
        % Calculate middle band moving average using convolution.
        cmid    = conv(nnandata, wtsvec);
        nnanmid = cmid(period:length(nnandata));
        
        % Calculate shift for the upper and lower bands. The shift is a
        % moving standard deviation of the data.
        mstd = nnandata(period:end); % Pre-allocate
        for i1 = period:length(nnandata)
           mstd(i1-period+1, :) = std(nnandata(i1-period+1:i1));
        end
        
        % Calculate the upper and lower bands.
        nnanuppr = nnanmid + nstd.*mstd;
        nnanlowr = nnanmid - nstd.*mstd;
        
        % Return the values.
        nanVec = nan(period-1,1);
        mid(~isnan(cl))  = [nanVec; nnanmid];
        uppr(~isnan(cl)) = [nanVec; nnanuppr];
        lowr(~isnan(cl)) = [nanVec; nnanlowr];
        
        % Format output
        vout = [mid,uppr,lowr];
        
    case 'keltner'  % Keltner Channels
        % Input data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable argument input
        if isempty(varargin)
            emaper = 20;
            atrmul = 2;
            atrper = 10;
        else
            emaper = varargin{1};
            atrmul = varargin{2};
            atrper = varargin{3};
        end
        
        % True range
        h_m_l = hi-lo;                                   % high - low
        h_m_c = [0;abs(hi(2:observ)-cl(1:observ-1))];  % abs(high - close)
        l_m_c = [0;abs(lo(2:observ)-cl(1:observ-1))];   % abs(low - close)
        tr = max([h_m_l,h_m_c,l_m_c],[],2);                 % true range
        
        % Average true range
        atr = ema(tr,atrper,observ);
        
        % Middle/Upper/Lower bands of keltner channels
        midd = ema(cl,emaper,observ);
        uppr = midd+atrmul*atr;
        lowr = midd-atrmul*atr;
        
        % Format output
        vout = [midd,uppr,lowr];
        
    case 'atr'      % Average True Range
        % Input data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable argument input
        if isempty(varargin)
            period = 20;
        else
            period = varargin{1};
        end
        
        % True range
        h_m_l = hi-lo;                                   % high - low
        h_m_c = [0;abs(hi(2:observ)-cl(1:observ-1))];  % abs(high - close)
        l_m_c = [0;abs(lo(2:observ)-cl(1:observ-1))];   % abs(low - close)
        tr = max([h_m_l,h_m_c,l_m_c],[],2);                 % true range
        
        % Average true range
        atr = ema(tr,period,observ);
        
        % Format Output
        vout = atr;
        
    case 'vr'       % Volatility Ratio
        % Input data
        hi = vin(:,1);
        lo = vin(:,2);
        cl = vin(:,3);
        
        % Variable argument input
        if isempty(varargin)
            period = 14;
        else
            period = varargin{1};
        end
        
        % True range
        h_m_l = hi-lo;                                   % high - low
        h_m_c = [0;abs(hi(2:observ)-cl(1:observ-1))];  % abs(high - close)
        l_m_c = [0;abs(lo(2:observ)-cl(1:observ-1))];   % abs(low - close)
        tr = max([h_m_l,h_m_c,l_m_c],[],2);                 % true range
        
        % Volatility Ratio
        vr = tr./ema(tr,period,observ);
        
        % Format Output
        vout = vr;
        
    case 'hhll'     % Highest High, Lowest Low
        % Input data
        hi = vin(:,1);
        lo = vin(:,2);
        
        % Variable argument input
        if isempty(varargin)
            period = 20;
        else
            period = varargin{1};
        end
        
        % Lowest Low
        llv = nan(observ,1);                        % preallocate lowest low
        llv(1:period) = min(lo(1:period));         % lowest low of first kperiods
        for i1 = period:observ                      % cycle through rest of data
            llv(i1) = min(lo(i1-period+1:i1));     % lowest low of previous kperiods
        end
        
        % Highest High
        hhv = nan(observ,1);                        % preallocate highest high
        hhv(1:period) = max(hi(1:period));        % highest high of first kperiods
        for i1 = period:observ                      % cycle through rest of data
            hhv(i1) = max(hi(i1-period+1:i1));    % highest high of previous kperiods
        end
        
        % Midpoint
        mp = (hhv+llv)/2;
        
        % Format Output
        vout = [hhv llv mp];
        
%--------------------------------------------------------------------------
        
    %%% Other
%==========================================================================
    case 'zigzag'   % ZigZag
        % Input data
        cl = vin;
        
        % Variable argument input
        if isempty(varargin)
            moveper = 7;
        else
            moveper = varargin{1};
        end
        
        % Preallocate zigzag
        zigzag = nan(observ,1);
        
        % First zigzag is first data point of input vector
        zigzag(1) = 1;
        
        i1 = 1; % index of input data
        i2 = 1; % index of output vector
        
        % The number of outputs is unknown
        while 1
            % Find the first value in the input, from the current index to
            % the number of observations, that has a price movement of
            % moveper greater or less than the value of current index of
            % the input and return the index of that value
            % If all of the following conditions are met, temp_var1 is the
            % index a zigzag
            temp_var1 = find(cl(i1:observ) > cl(i1)+cl(i1)*moveper/100 | ...
                cl(i1:observ) < cl(i1)-cl(i1)*moveper/100,1,'first');
            
            % If no value is found
            if     isempty(temp_var1)
                % If the current index is less than the number of
                % observations
                if i1 < observ
                    % If the current index of the output vector is greater
                    % than 1
                    if i2 > 1
                        % If the value of the input of the last recorded
                        % index is less than the value of the input of the
                        % index 2 recordings ago and there is a value
                        % between the value of the last recorded index and
                        % the number of observations that is less than the
                        % value of the last recorded index
                        if     cl(zigzag(i2)) < cl(zigzag(i2-1)) && ...
                                min(cl(i1:observ)) < cl(zigzag(i2))
                            % Find the index of the minimum value that is
                            % between the index of the last recorded value
                            % and the number of observations
                            [~,temp_var1] = min(cl(zigzag(i2):observ));
                            
                            % Set the output of the current index equal to
                            % the previously calculated index
                            zigzag(i2) = i1+temp_var1-1;
                            
                        % The opposite of the previous if statement
                        elseif cl(zigzag(i2)) > cl(zigzag(i2-1)) && ...
                                max(cl(i1:observ)) > cl(zigzag(i2))
                            [~,temp_var1] = max(cl(zigzag(i2):observ));
                            zigzag(i2) = i1+temp_var1-1;
                            
                        % The previous 2 statements are not true
                        % The output vector is complete
                        else
                            break
                        end
                        
                    % The previous statement is not true
                    % The output vector is complete
                    else
                        break
                    end
                    
                % The previous statement is not true
                % The output vector is complete
                else
                    break
                end
                
            % If the current index of the output vector is greater than 1
            elseif i2 > 1
                % If the value of the index of temp_var1 is greater than
                % the value of the last recorded index and the the value of
                % the last recorded index is greater than the value of the
                % index 2 recordings ago
                if     cl(temp_var1+i1-1) > cl(zigzag(i2)) && ...
                        cl(zigzag(i2)) > cl(zigzag(i2-1))
                    % Set the output of the current index equal to the
                    % index temp_var1
                    zigzag(i2) = temp_var1+i1-1;
                    
                % The opposit of the previous if statement
                elseif cl(temp_var1+i1-1) < cl(zigzag(i2)) && ...
                        cl(zigzag(i2)) < cl(zigzag(i2-1))
                    zigzag(i2) = temp_var1+i1-1;
                    
                % If the value of the input of the last recorded index is
                % less than the value of the input of the index 2 
                % recordings ago and there is a value between the value of
                % the last recorded index and temp_var1 that is less than
                % the value of the last recorded index
                elseif cl(zigzag(i2)) < cl(zigzag(i2-1)) && ...
                           min(cl(zigzag(i2):temp_var1+i1-1)) < cl(zigzag(i2))
                    % Find the index of the minimum value that is between
                    % the index of the last recorded value and temp_var1
                    [~,temp_var1] = min(cl(zigzag(i2):temp_var1+i1-1));
                    
                    % Set the output of the current index equal to the
                    % previously calculated index
                    zigzag(i2) = temp_var1+i1-1;
                    
                % The opposite of the previous statement
                elseif cl(zigzag(i2)) > cl(zigzag(i2-1)) && ...
                           max(cl(zigzag(i2):temp_var1+i1-1)) > cl(zigzag(i2))
                    [~,temp_var1] = max(cl(zigzag(i2):temp_var1+i1-1));
                    zigzag(i2) = temp_var1+i1-1;
                    
                % The previous 4 statements are not true
                else
                    % Increment the index of the output vector
                    % set the output of the incremented index equal to
                    % temp_var1
                    i2 = i2+1;
                    zigzag(i2) = temp_var1+i1-1;
                end
                
            % The current index of the output is equal to 1
            else
                % increment the index of the output vector
                % set the output of the incremented index equal to
                % temp_var1
                i2 = i2+1;
                zigzag(i2) = temp_var1+i1-1;
            end
            
            % Increment the index of the input data
            i1 = temp_var1+i1-1;
        end
        
        % Redefine the output data equal to the index of each zigzag, and
        % the value of the index of each zigzag
        zigzag = [zigzag(~isnan(zigzag)),cl(zigzag(~isnan(zigzag)))];
        
        % Format output
        vout = zigzag;
        
    case 'compare'  % Price Comparison
        % Input data
        numvars = size(vin,2);
        price   = vin;
        
        % Percent change relative to first price
        delta_percent = nan(observ,numvars);
        delta_percent(2:observ,:) = 100*(price(2:observ,:)-price(1,:))./price(1,:);
        
        % Format output
        vout = delta_percent;
        
    case 'pivot'        % Pivot Points
        % Input Data
        dt = vin(:,1);
        op = vin(:,2);
        hi = vin(:,3);
        lo = vin(:,4);
        cl = vin(:,5);
        
        % Variable Argument Input
        if isempty(varargin)
            type = 's';
        else
            type = varargin{1};
        end
        
        % Convert Matlab time to years, months, and days
        [year,month,day,~,~,~] = datevecmx(dt);
        
        % Frequency
        freq = diff(dt);
        if     sum(freq)/observ < 1 % Intraday
            freq = day;
        elseif sum(freq)/observ < 7 % Daily
            freq = month;
        else                        % Weekly/Monthly
            freq = year;
        end
        
        % Reassign open, high, low, and close based on frequency
        temp_var1 = unique(freq);
        num_dates = length(temp_var1);
        new_open  = nan(observ,1);
        new_high  = nan(observ,1);
        new_low   = nan(observ,1);
        new_close = nan(observ,1);
        for i1 = 2:num_dates
            last_per = freq == temp_var1(i1-1);
            this_per = freq == temp_var1(i1);
            temp_var2 = op(last_per);
            new_open (this_per) = temp_var2(1);
            new_high (this_per) = max(hi(last_per));
            new_low  (this_per) = min(lo(last_per));
            temp_var2 = cl(last_per);
            new_close(this_per) = temp_var2(end);
        end
        
        % Pivot Point
        switch type
            case 's'    % Standard
                pivot     = (new_high+new_low+new_close)/3;
                sprt(:,1) = pivot*2-new_high;
                sprt(:,2) = pivot-(new_high-new_low);
                res(:,1)  = pivot*2-new_low;
                res(:,2)  = pivot+new_high-new_low;
            case 'f'    % Fibonacci
                pivot     = (new_high+new_low+new_close)/3;
                sprt(:,1) = pivot-0.382*(new_high-new_low);
                sprt(:,2) = pivot-0.612*(new_high-new_low);
                sprt(:,3) = pivot-(new_high-new_low);
                res(:,1) = pivot+0.382*(new_high-new_low);
                res(:,2) = pivot+0.612*(new_high-new_low);
                res(:,3) = pivot+(new_high-new_low);
            case 'd'    % Demark
                X = nan(observ,1);
                temp_var1 = new_high+2*new_low+new_close;
                temp_var2 = 2*new_high+new_low+new_close;
                temp_var3 = new_high+new_low+2*new_close;
                X(new_close < new_open)  = temp_var1(new_close < new_open);
                X(new_close > new_open)  = temp_var2(new_close > new_open);
                X(new_close == new_open) = temp_var3(new_close == new_open);
                pivot = X/4;
                sprt  = X/2-new_high;
                res   = X/2-new_low;
        end
        
        % Format Ouput
        vout = [pivot sprt res];
        
    case 'sar'      % Parabolic SAR (Stop And Reverse)
        % Input Data
        hi = vin(:,1);
        lo = vin(:,2);
        
        % Variable Argument Input
        if isempty(varargin)
            step    = 0.02;
            maximum = 0.2;
        else
            step    = varargin{1};
            maximum = varargin{2};
        end
        af = step;
        
        % Directional Movement
        h_m_h = hi(2)-hi(1);                    % high - high
        l_m_l = lo(1)-lo(2);                      % low - low
        pdm   = 0;                                  % preallocate pdm1
        mdm   = 0;                                  % preallocate mdm1
        max_h = max(h_m_h,0);                       % max high
        max_l = max(l_m_l,0);                       % max low
        pdm(h_m_h > l_m_l) = max_h(h_m_h > l_m_l);  % +DM
        mdm(l_m_l > h_m_h) = max_l(l_m_l > h_m_h);  % -DM
        
        % false is long true is short
        new_dir            = false;
        new_dir(mdm < pdm) = true;
        
        % Defaults
        out_sar       = nan(observ,1);
        ep (new_dir)  = hi(2);
        sar(new_dir)  = lo(1);
        ep (~new_dir) = lo(2);
        sar(~new_dir) = hi(1);
        
        for i1 = 1:observ-1
            if new_dir
                % Switch to short if the low penetrates the SAR value
                if lo(i1+1) <= sar
                    new_dir = false;
                    sar = ep;
                    
%                     sar(sar < high(i1))   = high(i1);
%                     sar(sar < high(i1+1)) = high(i1+1);
                    
                    out_sar(i1+1) = sar;
                    
                    af = step;
                    ep = lo(i1+1);
                    
                    sar = sar+af*(ep-sar);
                    
%                     sar(sar < high(i1))   = high(i1);
%                     sar(sar < high(i1+1)) = high(i1+1);
                else
                    out_sar(i1+1) = sar;
                    
                    af(hi(i1+1) > ep) = af+step;
                    ep(hi(i1+1) > ep) = hi(i1+1);
                    af(af > maximum)    = maximum;
                    
                    sar = sar+af*(ep-sar);
                    
%                     sar(sar > low(i1))   = low(i1);
%                     sar(sar > low(i1+1)) = low(i1+1);
                end
            else
                % Switch to long if the high penetrates the SAR value
                if hi(i1+1) >= sar
                    new_dir = true;
                    sar = ep;
                    
%                     sar(sar > low(i1))   = low(i1);
%                     sar(sar > low(i1+1)) = low(i1+1);
                    
                    out_sar(i1+1) = sar;
                    
                    af = step;
                    ep = hi(i1+1);
                    
                    sar = sar+af*(ep-sar);
                    
%                     sar(sar > low(i1))   = low(i1);
%                     sar(sar > low(i1+1)) = low(i1+1);
                else
                    out_sar(i1+1) = sar;
                    
                    af(lo(i1+1) < ep) = af+step;
                    ep(lo(i1+1) < ep) = lo(i1+1);
                    af(af > maximum)   = maximum;
                    
                    sar = sar+af*(ep-sar);
                    
%                     sar(sar < high(i1))  = high(i1);
%                     sar(sar < high(i1+1)) = high(i1+1);
                end
            end
        end
        
        % Format Output
        vout = out_sar;
%--------------------------------------------------------------------------

end

end

%%% Simple Moving Average
%==========================================================================
function [vout status] = sma(vin,lag,observ)
% Set status
status = 1;

% If the lag is greater than or equal to the number of observations
if lag >= observ
    % End function, set status
    status = 0;
    return
end

% Preallocate a vector of nan's
vout = nan(observ,1);

% Simple moving average
ma = filter(ones(1,lag)/lag,1,vin);

% Fill in the nan's
vout(lag:end) = ma(lag:end);

end
%--------------------------------------------------------------------------

%%% Exponential Moving Average
%==========================================================================
function [vout status] = ema(vin,lag,observ)

% Preallocate output
vout   = nan(observ,1);

% Set status
status = 1;

% If the lag is greater than or equal to the number of observations
if lag >= observ
    status = 0;
    return
end

% Calculate the exponential percentage
k = 2/(lag+1);

% Calculate the simple moving average for the first 'exp mov avg' value.
vout(lag) = sum(vin(1:lag))/lag;

% K*vin; 1-k
kvin = vin(lag:observ)*k;
oneK = 1-k;

% First period calculation
vout(lag) = kvin(1)+(vout(lag)*oneK);

% Remaining periods calculation
for i1 = lag+1:observ
    vout(i1) = kvin(i1-lag+1)+(vout(i1-1)*oneK);
end

end
%--------------------------------------------------------------------------
