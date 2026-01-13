function result = approximation_to_GLd_eigenvalue (Z,d,n)
Id = eye(d);

one_matrix = ones(d);
A = ones(n) - eye(n);
ZA = Z .* kron(A,one_matrix);
deg = sum(A , 2);
D = diag(deg);

M = kron(D,Id)\ZA;
[V, L] = eig(M); 
[~,idx] = sort(real(diag(L)), 'descend');
result = V(:,idx(1:d));

end
