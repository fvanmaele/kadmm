%% KADMM for TV restoration
% Solves the image deblurring problem with spatially variant blur in Krylov
% subspace, using total variation as regularizer and ADMM as iterative method. 
%
% $$\min_{u\in
% K^{m}}TV(u)+\delta_{B(\varepsilon,f)}(Au)+\delta_{[0,1]^{n}}(u)$$
%
% The Krylov subspace is constructed using the Arnoldi method,
% with initial value set to $u0 = 0$. The vlaue of $\epsilon$ can be set
% with Morozov's discrepancy principle, $\epsilon = \sqrt{\tau MN}\sigma$,
% where a rule of thumb is $\tau = 1$, and \sigma^2 is the noise variance.
%
%% Input
% A:            forward operator for spatially variant blur
% f:            observed image, compressed into column vector
% V:            initial values for constrained variables (MNx3 matrix)
% D:            initial values for dual variables (MNx3 matrix)
% m:            order of Krylov subspace K(A, f)
% mu:           regularization parameter for total variation problem
% radius:       threshold for the residual ||Af - u||
% threshold:    threshold for relative change in each step
% max_iter:     maximum amount of ADMM iterations

%% Output
% u:            deblurred image, compressed into column vector

%% References
% File:         kadmmTV.m
% Author:       Ferdinand Vanmaele, University of Heidelberg
% Project:      TV restoration of spatially variant blur
%%
function [u, iters] = kadmmTV(A, f, V, D, m, mu, radius, threshold, max_iter)
    verbose = 1;

    % Matrix checks
    if ~ismatrix(A)
        error("A must be a matrix")
    end
    MN = size(A, 1);
    if MN ~= size(A, 2)
        error("A must be a quadratic matrix")
    end
    if MN ~= size(V, 1) || size(V, 2) ~= 3
        error("V must be an MNx3 matrix")
    end
    if MN ~= size(D, 1) || size(D, 2) ~= 3
        error("D must be an MNx3 matrix")
    end

    % Vector checks
    if ~isvector(f)
        error("f must be a vector")
    end
    if length(f) ~= MN
        error(['Sizes of A and f are not compatible, A = ' ...
                MN 'x' MN 'f = ' length(f) 'x1'])
    end
    
    % Scalar checks
    if ~isscalar(m) || m <= 0 || floor(m) ~= m
        error("m must be a positive integer")
    end
    if ~isscalar(mu)
        error("Mu must be a scalar")
    end   
    if ~isscalar(threshold) || threshold <= 0
        error("threshold must be a positive scalar")
    end
    if ~isscalar(max_iter) || max_iter <= 0
        error("max_iter must be a positive scalar")
    end

    % Orthonormal basis of m-th Krylov subspace
    u0 = zeros(length(f), 1);
    [Q, H] = arnoldi(A, f, m, u0);
    Qm = Q(:, 1:end-1);
    
    % Initial values for constrained variable v = Au
    v1 = V(:, 1);
    v2 = V(:, 2);
    v3 = V(:, 3);
 
    % Initial values for dual variable d
    d1 = D(:, 1);
    d2 = D(:, 2);
    d3 = D(:, 3);
    
    % Perform ADMM iteration
    I = eye(m);
    a = zeros(m, 1);
    k = 1;
    QHa = zeros(MN, 1);

    while k <= max_iter
        QHa_prev = QHa;
        % Solve least squares problem
        q = vertcat(Q'  * (v1 + d1), ...
                    Qm' * (v2 + d2), ...
                    Qm' * (v3 + d3));
        a = vertcat(H, I, I) \ q;

        % Matrix-vector products for ADMM
        QHa = Q * H * a; % O(n)
        Qa  = Qm * a;

        % Use relative change for termination
        if k >= 2
            relres = norm(QHa - QHa_prev) / norm(QHa - f);
            if mod(k, 250) == 0
                if verbose > 0
                    fprintf('relchange: %2.10f, iter: %d\n', relres, k)
                end
            end
            if relres <= threshold
                if verbose > 0
                    fprintf('relchange: %2.10f, iter: %d\n', relres, k)
                end
                break
            end
        end
        
        % Compute proximal operators
        v1 = projectToBall(QHa - d1, f(:), radius, 2);
        v2 = TV1D_denoise_tautString_mex(Qa - d2, mu);
        v3 = projectToBox(Qa - d3, 0, 1);
        
        % Update dual variables
        d1 = d1 - (QHa - v1);       
        d2 = d2 - (Qa - v2);
        d3 = d3 - (Qa - v3);

        % Proceed to next iteration
        k = k+1;
    end
    u = Qm * a;
    iters = k;
end

%%
function y = projectToBall(s, c, radius, type)
    nm = norm(s - c, type);
    
    if nm > radius
        y = c + radius * (s-c) / nm;
    else
        y = s;
    end
end

%%
function y = projectToBox(s, lower, upper)
    box_l = repelem(lower, length(s));
    box_u = repelem(upper, length(s));
    
    y = max(s, box_l');
    y = min(y, box_u');
end
