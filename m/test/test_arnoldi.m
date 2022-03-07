%% Test Arnoldi-MGS
N = 6;
p = 3;
A = unifrnd(0, 1, N, N);
b = unifrnd(0, 1, N, 1);
x0 = zeros(6, 1);
[Q, H] = arnoldi(A, b, p, x0);

%% Dimension checks
assert(size(Q, 1) == N);
assert(size(Q, 2) == p+1);
assert(size(H, 1) == p+1);
assert(size(H, 2) == p);

%% Orthonormality condition
[xx, yy] = meshgrid(1:p, 1:p);
idx = [xx(:) yy(:)];
eps = 1e-14;
dp = zeros(p, p);

for k = 1:size(idx, 1)
    [i, j] = deal(idx(k, 1), idx(k, 2));
    dp(i, j) = dot(Q(:, i), Q(:, j));
end
assert(norm(dp-eye(p)) <= eps)
