clear; clc

%% Blur and noise parameters
L = 21;
gamma = 4;
boundary = 'periodic';
BSNR = 60;

%% KADMM parameters
max_iter = 10000;
mu = 1;
order = 16;
TOL = 10e-6;
radf = @(n, var) 0.8544 * sqrt(n * var); % radius function

%% Evaluate performance
%u_orig = loadData(1);
u_orig = im2double(imread('001.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'corners', BSNR);

%%
fprintf("001, corners\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '001.jpg', 'JPEG');
imwrite(u_dmg, '001_60db_corners_sv.jpg', 'JPEG');
imwrite(u, '001_60db_corners_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(1);
u_orig = im2double(imread('001.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'center', BSNR);

%%
fprintf("001, center\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '001.jpg', 'JPEG');
imwrite(u_dmg, '001_60db_center_sv.jpg', 'JPEG');
imwrite(u, '001_60db_center_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(1);
u_orig = im2double(imread('001.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'gradient', BSNR);

%%
fprintf("001, gradient\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '001.jpg', 'JPEG');
imwrite(u_dmg, '001_60db_gradient_sv.jpg', 'JPEG');
imwrite(u, '001_60db_gradient_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(50);
u_orig = im2double(imread('050.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'corners', BSNR);

%%
%filename = sprintf("%3d_%ddb_%s", 50, BSNR, 'corners');
fprintf("050, corners\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '050.jpg', 'JPEG');
imwrite(u_dmg, '050_60db_corners_sv.jpg', 'JPEG');
imwrite(u, '050_60db_corners_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(50);
u_orig = im2double(imread('050.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'center', BSNR);

%%
fprintf("050, center\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '050.jpg', 'JPEG');
imwrite(u_dmg, '050_60db_center_sv.jpg', 'JPEG');
imwrite(u, '050_60db_center_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(50);
u_orig = im2double(imread('050.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'gradient', BSNR);

%%
fprintf("050, gradient\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '050.jpg', 'JPEG');
imwrite(u_dmg, '050_60db_gradient_sv.jpg', 'JPEG');
imwrite(u, '050_60db_gradient_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(103);
u_orig = im2double(imread('103.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'corners', BSNR);

%%
fprintf("103, corners\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '103.jpg', 'JPEG');
imwrite(u_dmg, '103_60db_corners_sv.jpg', 'JPEG');
imwrite(u, '103_60db_corners_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(103);
u_orig = im2double(imread('103.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'center', BSNR);

%%
fprintf("103, center\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '103.jpg', 'JPEG');
imwrite(u_dmg, '103_60db_center_sv.jpg', 'JPEG');
imwrite(u, '103_60db_center_sv_restored.jpg', 'JPEG');

%%
%u_orig = loadData(103);
u_orig = im2double(imread('103.jpg'));
[MSE, u, u_dmg, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, 'gradient', BSNR);

%%
fprintf("103, gradient\n");
fprintf("MSE: %f, radius: %f, iterations: %d, time: %f\n", MSE, radius, iters, elapsed);
%imwrite(u_orig, '103.jpg', 'JPEG');
imwrite(u_dmg, '103_60db_gradient_sv.jpg', 'JPEG');
imwrite(u, '103_60db_gradient_sv_restored.jpg', 'JPEG');
