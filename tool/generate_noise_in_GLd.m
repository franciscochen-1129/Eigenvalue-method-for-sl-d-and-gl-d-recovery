function Z = generate_noise_in_GLd(G,d,n,theta)
Z = eye(d*n);

for i = 1:n
    rows_i = (i-1)*d + 1 : i*d;
    G_i = G(rows_i, :);

    for j = 1:n
        rows_j = (j-1)*d + 1 : j*d;

        if j > i
            G_j = G(rows_j, :);

            rel = G_i \ G_j;
            X = randn(d);
            X = X / norm(X,'fro');

            Z(rows_i, rows_j) = rel * expm(theta * X);

        elseif j < i
            Z(rows_i, rows_j) = Z(rows_j, rows_i) \ eye(d);

        else
            continue;
        end
    end
end
end
