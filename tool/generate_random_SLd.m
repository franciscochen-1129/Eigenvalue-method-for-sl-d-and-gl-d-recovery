function G = generate_random_SLd(n,d,theta)
G = zeros(n*d,d);
Id = eye(d);
for t = 1:n
    row = (t - 1)*d + 1 : t*d;
    X = randn(d);
    X = X - (trace(X)/d) * Id;
    G(row, :) = expm(theta*X);
end