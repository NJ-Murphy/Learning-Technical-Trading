function b = zOnlineAnticor(Data, T, n1, H0)
% 
% n1: look-back period for the indicator
% T: index of the most recent data we have
% Data: returns/price relatives of the stocks over time (time x
% stocks)
% H0 are the controls for the previous period

    % Check if there is enough data
    if T < 2*n1
        b = H0; %if not return previous periods portfolio 
    else
        H0(~any(~isnan(H0), 2),:) = 0;

        N = length(H0);
        
        % Get the data for the 2 windows
        Window1 = Data((end-2*n1+1):(end-n1),:)-1; 
        Window2 = Data((end -n1+1):end,:)-1; 
        
        % Compute means for 2 windows
        mu1 = mean(Window1);
        mu2 = mean(Window2);
        claims = zeros(N);
        H1 = zeros(size(H0));

        if n1 > 1
            sigma1 = transpose(std(Window1,0,1));
            sigma2 = transpose(std(Window2,0,1));
            M_cov = (1/(n1 - 1))*(transpose(Window1-(ones(n1,1)*mu1)))*(Window2-(ones(n1,1)*mu2));
            M_cor = M_cov ./ (sigma1 *transpose(sigma2));
            M_cor(find(isnan(M_cor))) = 0;

            for j = 1:N
              claims(:,j) = Claim(j, M_cor(:,j),diag(M_cor), N,claims(:,j),mu2);
            end
        else
            sigma1 = zeros(N,1);
            sigma2 = zeros(N,1);
            M_cov = zeros(N);
            M_cor = ones(N);
            claims = eye(N);
        end

        for i = 1:N
          H1(i) = H0(i) + 1/3*(sum(claims(:,i) - claims(i,:)'));
        end

        b = H1';
        b = ReNormAgMix(b);

    end
end

function Claimj = Claim(j, Mcorj , McorDiag , M, Claimj , mu)
 for i = 1:M
    if mu(i) >= mu(j) && Mcorj(i) > 0
      Claimj(i) = Mcorj(i) + max (0, -McorDiag(i))+max(0, -McorDiag(j));
    end
 end
end

function q1 = ReNormAgMix(q0)
  absSum = sum(abs(q0 - mean(q0 )));
  if absSum > eps ()
    q1 = (q0 - mean (q0 ))/absSum;
  else
    q1 = zeros(size(q0));
  end
  %a = sum(abs .(q1));
end