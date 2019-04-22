function b = OnlineZBCRP (Data, n1)
    exp_tup = Data;
    M = size(exp_tup,2);
    onesM = ones(M,1);

    mu = transpose(mean(exp_tup ,1));
    
    if n1 < 3
      sigma = eye(M);
      inv_sigma = sigma;
    else
      sigma = cov(exp_tup);
      if cond(sigma) < (1/(10^(-2))) && n1 > M
        inv_sigma = inv(sigma);
      else
        sigma = diag(diag(sigma));
        if cond(sigma) > 1/(10^(-2))
          sigma = eye(M);
          inv_sigma = inv(sigma);
        elseif any(sigma== 0) == true
          sigma = ZeroVar(sigma);
          inv_sigma = inv(sigma);
        end
      end
    end
    sigma_scale = transpose(onesM)*inv_sigma*onesM;
    risk_aver = onesM'*abs((inv_sigma*(mu*transpose(onesM) - onesM *...
        transpose(mu))*inv_sigma*onesM )/sigma_scale);
    if risk_aver ==0
      risk_aver = 1;
    end
    b = inv_sigma/risk_aver*((mu*transpose(onesM) - onesM * transpose (mu ))/...
    ( sigma_scale ))*inv_sigma* onesM;
end

function S = ZeroVar(S)
  Z_Ind = ( diag(S) == 0);
  S(Z_Ind , Z_Ind) = mean(diag(S))*eye(sum(Z_Ind));
end