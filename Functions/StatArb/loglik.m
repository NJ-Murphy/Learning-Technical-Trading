function ll = loglik(x,deltanu)
  i = 1:length(deltanu);
  % mu = x(1), lambda = x(2), sigma2 = x(3)
  ll = (1/2)*sum(log(x(3)*(i.^(2*x(2))))) + (1/(2*x(3)))*sum(i.^(-2*x(2))*(deltanu - x(1)).^2);
  
end