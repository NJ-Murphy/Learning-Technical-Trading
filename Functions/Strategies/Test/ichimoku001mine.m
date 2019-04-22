function out = ichimoku(Data, n1, n2)
%varargin)
    % ICHIMOKU Ichimoku Kinko Hyo (Cloud Charts) indicators
    %  
    % ICHIMOKU(P) Plot the Ichimoku Kinko Hyo using prices P as an NxM
    %   matrix with N days and M stocks. The Ichimoku Kinko Hyo ("At a glance 
    %   equilibrium chart") consists of five lines and the Kumo (The Cloud). The 
    %   calculation for four of these lines involves taking only the midpoints 
    %   of previous highs and lows. Ichimoku uses three key time periods for 
    %   its input parameters: 7, 22, and 44. This was historically 9,26 and 52 
    %   as based on a 6 day week. P can be a FINTS object containing the same
    %   inputs as required for CANDLE.
    %
    %   1. Tenkan-Sen    (Conversion Line)
    %   2. Kijun-Sen     (Base Line)  
    %   3. Chikou Span   (Lagging Span)  
    %   4. Senkou Span A (Leading Span A)
    %   5. Senkou Span B (Leading Span B)
    % 
    % The five lines are calculated as follows:
    % +-----------------+----------------------------------+------+----------+
    % | Line            | Calculation                      | Per. | Lag/Lead |
    % +-----------------+----------------------------------+------+----------+
    % | Conversion Line | (Highest High + Lowest Low)/2    | 7    | -        |
    % | Base Line       | (Highest High + Lowest Low)/2    | 22   | -        |
    % | Lagging Span    | Closing Price (lagged)           | 1    | lag 22   |
    % | Leading Span 1  | (Conversion Line + Base Line )/2 | 1    | lead 22  |
    % | Leading Span 2  | (Highest High + Lowest Low)/2    | 44   | lead 22  |
    % | Cloud           | Between Leading Span 1 and 2     | -    | -        |
    % +-----------------+----------------------------------+------+----------+
    %
    % ICHIMOKU(OP,HI,LO,CL) Plot Ichimoku chart using the opening price (OP), 
    %   high price (HP), the low price (LP) and the closing price (CP). Price 
    %   is plotted as a candel stick.
    %
    % [L,I,T] = ICHIMOKU(P) Compute Ichimoku lines as NxMx5 matrix LINES 
    %   for N days and M stocks using prices P. The buy-sell indicator S 
    %   ranging between +4 and -4 from strong buy to strong sell as an NxM 
    %   matrix. The trend is computed +1 (-1) for up (down) as an NxM matrix.
    %
    % Buy and sell signals are given with the crossover technique
    % +---------+-----------------+------------------------------------------+
    % | Signal  | Cross-over      | Description                              |
    % +---------------------------+------------------------------------------+
    % | ++Buy   | Strong Bullish  | Bullish with price above the Cloud       |
    % |  +Buy   | Bullish         | Conversion crosses Base Line from below  |
    % |   Buy   | Normal Bullish  | Bullish cross-over within the Cloud      |
    % |  -Buy   | Weak Bullish    | Bullish with price below the Cloud       |
    % |  +Sell  | Weak Bearish    | Bearish with price above the Cloud       |
    % |   Sell  | Normal Bearish  | Bearish cross-over within the Cloud      |
    % |  -Sell  | Bearish         | Conversion crosses Base Lind from above  |
    % | --Sell  | Strong Bearish  | Bearish with price below the Cloud       | 
    % +---------+-----------------+------------------------------------------+
    %
    %   If the price is above (below) the Cloud the prevailing trend is said 
    %   to be up (down). If the Lagging Span was below the closing price and a 
    %   sell signal was issued, then the strength is with the sellers, 
    %   otherwise it is a weak signal. Support and resistance levels are
    %   provided by the Cloud. Price falling into the Cloud from above are 
    %   supported and prices climbing into the cloud from below are resisted. 
    %
    % [LINES,I,T] = ICHIMOKU(OP,HI,LO,CL) compute lines, indicators and trend
    %   using opening, high, low and closing prices.
    %
    % [LINES,I,T] = ICHIMOKU(OP,HI,LO,CL,PL) compute lines, indicators and 
    %   trend using opening, high, low and closing prices. PL is a 3x1 vector
    %   with the user defined time-periods. The default value is [7,22,44].
    %
    % Example 1:
    %   >> p = cumsum(randn(100,3)* 0.10);
    %   >> ichimoku(p);
    %
    % Example 2:
    %   >> p = cumsum(randn(100,3)* 0.10);
    %   >> [f,i,t]=ichimoku(p);
    % 
    % See Also: 

    % Authors: Tim Gebbie

    % FIXME: FINTS 4xDATES, FINTS N stocks all close prices, NxM matrix

    % $Id: ichimoku.m 112 2011-08-22 14:01:41Z tgebbie $

    %% Inputs
   % pl = [7,22,44];
    %switch nargin
    %    case 1
    %        if isa(varargin{1},'fints')
    %            fts = varargin{1};
    %            op = fts2mat(fts.Open);
    %            hi = fts2mat(fts.High);
    %            lo = fts2mat(fts.Low);
    %            cl = fts2mat(fts.Close);
    %        else
    %            op = varargin{1};
    %            hi = op;
    %            lo = op;
    %            cl = op;
    %        end
    %    case 4
    %            op = varargin{1};
    %            hi = varargin{2};
    %            lo = varargin{3};
    %            cl = varargin{4};
    %   case 5
    %            op = varargin{1};
    %            hi = varargin{2};
    %            lo = varargin{3};
    %            cl = varargin{4};
    %            pl = varargin{5};
    %    otherwise
    %        error('Incorrect Input Arguments');
    %end
    
    % Redefine the function inputs
    n1 = Data{2};
    n2 = Data{3};
    Data = reshape(Data{1},size(Data{1},2),size(Data{1},3))';
   % fts(:,1) = Data(:,2);
   % fts(:,2) = Data(:,3);
   % fts(:,3) = Data(:,3);
   % fts(:,4) = Data(:,4);
    
    op = Data(:,2);
    hi = Data(:,3);
    lo = Data(:,4);
    cl = Data(:,1);
    
    pl(1) = 7;
    pl(2) = n1;
    pl(3) = n2;
    %% Compute the 5 lines
    m = size(Data);
    f = NaN * zeros(m(1)+pl(2),5,m(2));
    for i=1:m(1)
        % 1. Conversion Line
        if i>pl(1)
            f(i,1,:) =  (max(hi(i-pl(1):i,:)) + min(lo(i-pl(1):i,:)))/2;
        end
        % 2. Base Line
        if i>pl(2)
            f(i,2,:) = (max(hi(i-pl(2):i,:)) + min(lo(i-pl(2):i,:)))/2;
        end
        % 3. Lagging Span
        if i>pl(2)
            f(i-pl(2),3,:) = cl(i,:);
        end
        % 4. Leading Span 1
        f(i+pl(2),4,:) = (f(i,1,:) + f(i,2,:))/2;
        % 5. Leading Span 2
        if i>pl(3) 
            f(i+pl(2),5,:) = (max(hi(i-pl(3):i,:)) + min(lo(i-pl(3):i,:)))/2;
        end  
    end
    
    s =  repmat(cl,1,5);
    %% Create the indicators
    % the bull indicator (+ if conversion is above base and crossed from below)
    bull = [NaN * zeros(1,5); squeeze(2*sign(diff(sign(f(1:m(1),1,:) - f(1:m(1),2,:)))))];
    %
    % cloud indicator
    for j=1:m(2)
        ac = s(:,j)>max(f(1:m(1),[4,5],j)')';
        bc = s(:,j)<min(f(1:m(1),[4,5],j)')';
        cloud(ac,j) = +1; % closing price above cloud
        cloud(bc,j) = -1; % closing price below cloud
    end
    % construct buy-sell indicators
    ind = cloud + bull;
    
    out{1} = ind(:,1);
    out{2} = f;
    out{3} = cloud;
end


