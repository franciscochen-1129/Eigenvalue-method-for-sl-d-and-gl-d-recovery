clear; clc;

addpath(fullfile(fileparts(mfilename('fullpath')),'tool'));

seed = 3;
rng(seed);

n = 10;
d = 3;
theta = 1; 
sigma = 1e-10; 

Id = eye(d);
one_matrix = ones(d);

G = generate_random_SLd(n,d,theta);
Z = generate_noise_in_GLd(G,d,n,sigma);
result = approximation_to_GLd_eigenvalue (Z,d,n);

Adj = ones(n) - eye(n);
deg = sum(Adj , 2);
D = diag(deg);
ZA = Z .* kron(Adj,one_matrix);

range = 1 : d;
T = result(range,:);
scale = det(T)^(1/d);
result = result./scale;


Dbig = sparse(kron(D, Id));
Lbig = Dbig - ZA;

error = norm(Lbig * result, 'fro')^2;

disp(error)