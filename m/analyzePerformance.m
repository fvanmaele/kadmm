%% Performance analysis of KADMM
% Creates a blurred and noisy version of an image, restores it using KADMM
% and computes the MSE between both images. All parameters except the
% initial values of the problem (u0, V, D) and the type of noise (Gaussian) 
% are set by the user. The forward operator only depends on the size of the
% given image, and additional blur parameters.
%
%% Input
% u_orig:       original image
% max_iter:     maximum amount of iterations for KADMM
% mu:           regularization parameter for (TV) proximal operator
% radf:         function which computes the radius ||Au-f|| < E
% order:        order of Krylov subspace
% TOL:          tolerance for relative change of residuals (termination)
% L:            size of the blurring kernel, LxL
% gamma:        parameter for the blur variance
% boundary:     type of boundary condition (periodic, zero, reflexive)
% blurType:     type of spatially variant blur (corners, center, gradient)
% BSNR:         blurred signal-to-noise ratio (if Inf, do not add noise)

%% Output
% u_restored:   restored image
% u_damaged:    image damaged by spatially variant blur and gaussian noise
% iters:        amount of ADMM iterations performed before terminating
% MSE:          mean-squared-error between original and restored image
% elapsed:      elapsed time in seconds

%% References
% File:         analyzePerformance.m
% Author:       Ferdinand Vanmaele, University of Heidelberg
% Project:      TV restoration of spatially variant blur
%%
function [MSE, u_restored, u_damaged, iters, elapsed, radius] = analyzePerformance(u_orig, max_iter, mu, radf, order, TOL, ...
    L, gamma, boundary, blurType, BSNR)

    % Precondition checks
    if ~ismatrix(u_orig)
        error("u_orig must be a matrix")
    end
    if ~isscalar(max_iter) || max_iter ~= floor(max_iter) || max_iter < 1
        error("max_iter must be a positive integer")
    end
    if ~isscalar(mu)
        error("mu must be a scalar")
    end
    if ~isscalar(order) || order ~= floor(order) || order < 1
        error("order must be a positive integer")
    end
    if ~isscalar(TOL) || TOL < 0
        error("TOL must be a positive scalar")
    end
    if ~isscalar(L) || L ~= floor(L) || mod(L, 2) ~= 1 || L < 1
        error("L must be a positive and uneven integer")
    end
    if ~isscalar(gamma)
        error("gamma must be a scalar")
    end
    if ~isscalar(BSNR) || BSNR < 0
        error("BSNR must be a positive scalar")
    end
    if ~ischar(boundary) || ~ischar(blurType)
        error("boundary and blurType must be character arrays")
    end
    if ~isa(radf, 'function_handle')
        error("radf must be a function handle")
    end
    
    % Compute image dimension
    height = size(u_orig, 1);
    width = size(u_orig, 2);
    pixels = height*width;
    
    % Compute blurred (spatially variant) and noisy image
    A = forward_Sv2d(height, width, L, gamma, boundary, blurType);
    B = reshape(A*u_orig(:), height, width);

    if BSNR == Inf
        u_damaged = B;
        radius = TOL;
    else
        noise_var = var(B(:)) / exp(BSNR / 10);
        u_damaged = imnoise(B, 'gaussian', 0, noise_var);
        radius = radf(pixels, noise_var);
    end
    V = unifrnd(0, 1, pixels, 3);
    D = unifrnd(0, 1, pixels, 3);
    
    % Run KADMM
    tic
    [u_restored, iters] = kadmmTV(A, u_damaged(:), V, D, order, mu, radius, TOL, max_iter);
    u_restored = reshape(u_restored, height, width);
    elapsed = toc;
    
    % Compute mean-square error
    MSE = immse(u_orig, u_restored);
end
