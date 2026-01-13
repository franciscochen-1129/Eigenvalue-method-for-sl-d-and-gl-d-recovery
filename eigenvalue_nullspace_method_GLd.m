clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')),'tool'));

seed = 3;
rng(seed);

n = 10;
d = 3;
theta = 2;
sigma = 0.01;

G = generate_random_GLd(n,d,theta);
Z = generate_noise_in_GLd(G,d,n,sigma);

Id = eye(d);
one_matrix = ones(d);

A = ones(n) - eye(n);
ZA = Z .* kron(A,one_matrix);

deg = sum(A , 2);
D = diag(deg);

M = kron(D,Id)\ZA;
[V, L] = eig(M);
[~,idx] = sort(real(diag(L)), 'descend');
V = V(:,idx);

result_eigenvalue_method = V(:,1:d);
result_eigenvalue_method = real(result_eigenvalue_method);

result_nullspace_method = null(kron(D,Id) - ZA);

Dbig = sparse(kron(D, Id));
Lbig = Dbig - ZA;

spec_error_eigen = norm(Lbig * result_eigenvalue_method, 'fro')^2;
spec_error_null  = norm(Lbig * result_nullspace_method, 'fro')^2;

% Using directly optimization
nd = n*d;

residue = @(X) norm((Dbig - ZA) * X,'fro')^2;
fun = @(x) residue(reshape(x, nd, d));

nonlcon = @(x) deal( [], reshape( ...
    (reshape(x, nd, d).'*reshape(x, nd, d) - n*Id), [], 1) );

X0 = real(result_eigenvalue_method);
[Q,~] = qr(X0,0);
X0 = sqrt(n)*Q;
x0 = X0(:);

minimizer = fmincon(fun,x0,[],[],[],[],[],[],nonlcon);

Xopt = reshape(minimizer, nd, d);
error_optimized = fun(minimizer);

disp(error_optimized);
disp(spec_error_eigen);
