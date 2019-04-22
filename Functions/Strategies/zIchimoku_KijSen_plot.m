function out = Ichimoku_KijSen_plot(Data, n1, n2)
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

    % Authors: Tim Gebbie, Nic Murphy

    % FIXME: FINTS 4xDATES, FINTS N stocks all close prices, NxM matrix

    % $Id: ichimoku.m 112 2011-08-22 14:01:41Z tgebbie $

    %% Inputs
    
    % Redefine the function inputs
    n1 = Data{2};
    n2 = Data{3};
    m = size(Data{1});
    Data = Data{1};
    
    
    op = Data(:,2,:);
    hi = Data(:,3,:);
    lo = Data(:,4,:);
    cl = Data(:,1,:);
    
    %this can be made faster by replacing below in above step
    op = reshape(op(:,1,:),size(op(:,1,:),1),size(op(:,1,:),3))';
    hi = reshape(hi(:,1,:),size(hi(:,1,:),1),size(hi(:,1,:),3))';
    lo = reshape(lo(:,1,:),size(lo(:,1,:),1),size(lo(:,1,:),3))';
    cl = reshape(cl(:,1,:),size(cl(:,1,:),1),size(cl(:,1,:),3))'; 
        
    pl(1) = 7;
    pl(2) = n1;
    pl(3) = n2;
    %% Compute the 5 lines
      %[m(1),~,m(2)] = size(Data);
    f = nan(5,m(3)+pl(2),m(1));
    for i = 1:m(3)
        % 1. Conversion Line/Tenkan-sen
        if i>pl(1)
            f(1,i,:) =  (max(hi(i-pl(1):i,:)) + min(lo(i-pl(1):i,:)))/2;
        end
        % 2. Base Line/Kijun-sen
        if i>pl(2)
            f(2,i,:) = (max(hi(i-pl(2):i,:)) + min(lo(i-pl(2):i,:)))/2;
        end
        % 3. Lagging Span/Chikou Span
        if i>pl(2)
            f(3,i-pl(2),:) = cl(i,:);
        end
        % 4. Leading Span 1/Senkou Span A
        f(4,i+pl(2),:) = (f(1,i,:) + f(2,i,:))/2;
       
        % 5. Leading Span 2/Senkou Span B
        if i>pl(3) 
            f(5,i+pl(2),:) = (max(hi(i-pl(3):i,:)) + min(lo(i-pl(3):i,:)))/2;
        end  
    end
    
    % Kijun Sen Cross
    sig = nan(size(f,2),1);
   for j = 1:m(1) 
    for i = (pl(2)+pl(3)+1):m(3) %time loop
        if (f(2,i,:)<=cl(i,j)) & (cl(i-1,j)<f(2,i-1,:))%bullish cross
            if (cl(i,j)>f(4,i,:)) & (cl(i,j)>f(5,i,:))
              sig(i,j) = +3;
            elseif f(4,i,:)<=cl(i,j) & cl(i,j)<=f(5,i,:)
              sig(i,j) = +2;
            elseif cl(i,j)< f(4,i,:) & cl(i,j)<f(5,i,:)
                sig(i,j) = +1;
            end
        elseif (f(2,i,:)>cl(i,j)) & (cl(i-1,j)>=f(2,i-1,:)) %bearish cross
            if (cl(i,j)<f(4,i,:)) & (cl(i,j)<f(5,i,:))
              sig(i,j) = -3;
            elseif f(4,i,:)<=cl(i,j) & cl(i,j)<=f(5,i,:)
              sig(i,j) = -2;
            elseif cl(i,j)> f(4,i,:) & cl(i,j)>f(5,i,:)
                sig(i,j) = -1;
            end
        else
            sig(i,j) = 0;
        end
    end
   end
    
    % convert signals into 1, -1 and 0 rather than strong, weak etc for
    % temporary use
%     sig = sig(m(3),:);
%     sig(sig>0) = 1;
%     sig(sig<0) = -1;
%     out{1} = sig;
    %% for plotting
    out{1} = sig(:,1);
    out{2} = f;
    out{3} = cl;
end