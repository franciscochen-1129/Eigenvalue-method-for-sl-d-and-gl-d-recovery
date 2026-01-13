function G = generate_random_GLd(n,d,theta)
G = zeros(n*d,d);
for t = 1:n
    row = (t - 1)*d + 1 : t*d;
    X = randn(d);
    G(row, :) = expm(theta*X);
end


