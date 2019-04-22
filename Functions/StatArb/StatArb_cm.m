function [Min_t,t_c_95,p_val_95,p_loss,Min_tsim] = StatArb_cm(S)
 
%compute discounted log wealth fluctuations
deltanu = diff(S);
n = length(deltanu);
i = 1:n;

%% Min-t for overall strategy
% minimize negative ll
lb = [-Inf,-Inf,0];
ub = [];
x0 = [0.1,0.1,0.1];
options = optimoptions('fmincon','Display','off');

[x,fval,~,~] = fmincon(@(x)loglik(x,deltanu),x0,[],[],[],[],lb,ub,[],options);

%compute the standard errors of the maximium likelihood procedure using optimHess
hess_analytic = analytical_hess(x(1),x(2),x(3),0,i,deltanu);

%standard errors
SE = sqrt(diag(inv(-hess_analytic)));

% compute the t-statistics
t_mu = x(1)/SE(1);
t_lam = -x(2)/SE(2);
%t_thetlamsim = (-x(1)+0.5)/SEsim(2);
%t_thetlam1sim = (-x(1))/SEsim(2);

Min_t = min(t_mu,t_lam);

%% Monte Carlo for critical values for CM model
%initialise parameters
numsim = 5000;
lamsim = 0;
musim = 0;
sigmasim = 0.01;

rng(5)

Min_tsim = zeros(0,numsim);
for m = 1:numsim % numsim simulations
    deltanusim = zeros(0,n);
    z = randn(n,1);
    for i = 1:n %simulate path of length n
      %simulate a single path for a given parameter combo and given innovation z
      deltanusim(i) = musim + sigmasim*z(i)*i^lamsim;
    end
    i = 1:n;
    [x,fval,~,~] = fmincon(@(x)loglik2(x,deltanusim),x0,[],[],[],[],lb,ub,[],options); %ML estimates of simulated path
    
    %compute the standard errors of the maximium likelihood procedure using optimHess
    hess_analytic = analytical_hess(x(1),x(2),x(3),0,i,deltanusim');

    %standard errors
    SEsim = sqrt(diag(inv(-hess_analytic)));
    
    % compute the 5 t-statistics
    t_musim = x(1)/SEsim(1);
    t_lamsim = -x(2)/SEsim(2);
    %t_thetlamsim = (-x(1)+0.5)/SEsim(2);
    %t_thetlam1sim = (-x(1))/SEsim(2);

    Min_tsim(m) = min(t_musim,t_lamsim);
    
end


%% Compute the critical t stats
t_c_95 = quantile(Min_tsim, 0.95);
t_c_90 = quantile(Min_tsim, 0.90);

t_c = Min_tsim(Min_tsim>t_c_95);
%t_c_95 <- max(Min_tsim)


%% compute p-value
p_val_95 = 1-sum(Min_t>Min_tsim(Min_tsim>t_c_95))/length(Min_tsim(Min_tsim>t_c_95));
p_val_90 = sum(Min_t > t_c_90)/length(Min_tsim);

%% Probability of loss per period

p_loss = zeros(length(n),1);
for k = 1:n
   i = 1:k;
  [x,fval,~,~] = fmincon(@(x)loglik(x,deltanu(1:k)),x0,[],[],[],[],lb,ub,[],options);

  p_loss(k) = normcdf((-x(1)*k*(k+1)/2)/(sqrt(x(3)*sum(i.^(2*x(2))))));
end


end