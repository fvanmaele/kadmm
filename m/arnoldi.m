%% Arnoldi-MGS
% Creates an orthonormal basis of the p-order Krylov subspace K_p(A, b-Ax0)
%
%% Input
% A:    coefficient matrix of the linear system Ax = b, assumed quadratic
% b:    right-hand side
% x0:   start approximation
% p:    order of the Krylov subspace
%
%% Output
% Q:    orthonormal basis of Krylov subspace
% H:    upper Hessenberg matrix such that A*Q[:,1:p] = Q*H
%
%% References
% File:       arnoldi.m
% Author:     Ferdinand Vanmaele, University of Heidelberg
% Project:    TV restoration of spatially variant blur
%%
function [Q, H] = arnoldi(A, b, p, x0)
    % Precondition checks
    if ~ismatrix(A)
        error("A must be a matrix")
    end

    if ~isvector(b)
        error("b must be a vector")
    end
    
    if ~isvector(x0)
        error("x0 must be a vector")
    end
    
    if length(b) ~= length(x0)
        error("b and x0 must have the same size")
    end
    
    if ~isscalar(p)
        error("p must be an integer")
    end
    
    [N, M] = size(A);
    if N ~= M
        error("A must be a quadratic matrix")
    end
    
    if length(b) ~= N
        error(['A and b have incompatible sizes, A = ' N 'x' m 'b = ' length(b) 'x1'])
    end
    
    H = zeros(p+1, p);
    Q = zeros(N, p);
    r0 = b - A*x0;
    Q(:, 1) = r0 / norm(r0);
    
    for i = 1:p
        v = A*Q(:, i);
        for j = 1:i
            H(j, i) = dot(Q(:, j), v);
            v = v - H(j, i)*Q(:, j);
        end
        H(i+1, i) = norm(v);
        Q(:, i+1) = v / H(i+1, i);
    end
end