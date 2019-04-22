function p = gen_port(n_portfolios, n_assets)
rng(10)
    function b =  generate_portfolio(n_assets)
        %Step 1:
        x = gamrnd(1,1,n_assets,1);
        %x = -1 + (2)*rand(n_assets,1);
        %Step 2:
        b = x/sum(x);
        b(n_assets) = 1-sum(b(1:(n_assets-1)));
    end

p = zeros(n_portfolios,n_assets);
for i = 1:n_portfolios
  p(i,:) = generate_portfolio(n_assets);
end
end