clear; clc;

%% Load data (spatially invariant)
[A, b, u_orig] = deblurringLoadData('BSDS500', 1, false);
imshow(u_orig);

%%
height = size(u_orig, 1);
width = size(u_orig, 2);
pixels = height*width;
B = reshape(b, height, width);
imshow(B);

%%
mu = 1;
order = 16;
radius = 1e-6;
threshold = 1e-5;
V = zeros(pixels, 3);
D = zeros(pixels, 3);

%%
[u, iters] = kadmmTV(A, b, V, D, order, mu, radius, threshold, 5000);
U = reshape(u, height, width);

%%
subplot(1,2,1); imshow(B); title('Degraded');
subplot(1,2,2); imshow(U); title('Restored');
