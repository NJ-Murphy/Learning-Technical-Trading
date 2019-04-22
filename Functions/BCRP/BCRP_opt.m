function [b0,fval] = BCRP_opt(n_assets,x)
% minimize negative ll
lb = zeros(n_assets,1);
ub = zeros(n_assets,1);
x0 = repmat(0.1,n_assets,1);
options = optimoptions('fmincon','Display','off','MaxIterations',5000);


[b0,fval,~,~] = fmincon(@(b)maxwealth(b,x),x0,[],[],ones(n_assets,1)',1,lb,ub,[],options);

end