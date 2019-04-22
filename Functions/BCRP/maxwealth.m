function S = maxwealth(b,x)
 jh = cumprod(x*b);
 S = -jh(end);
end