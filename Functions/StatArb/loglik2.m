% this log-likelihood is created for use in MC for different simulated
% incremtal processes
function ll2 = loglik2(x,deltanu)
  i = 1:length(deltanu);
  % mu = x(1), lambda = x(2),  sigma2 = x(3)
  %negative log-likelihood function
  ll2 = (1/2)*sum(log(x(3)*(i.^(2*x(2))))) + (1/(2*x(3)))*sum((1./(i.^(2*x(2))))*(deltanu' - x(1)).^2);
end