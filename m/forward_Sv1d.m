%% forward_Sv1d
% Forward operator for 2d separable Gaussian blur (1d component),
% discretized to size LxL
%
%% Input
% N:            size of 1d signal
% L:            kernel size
% sigma:        variance of gaussian blur, spatially variant (callable)
% boundary:     type of boundary condition
%
%% Output
% A:            forward operator
%
%% References
% File:         forward_Sv1d.m
% Author:       Ferdinand Vanmaele, University of Heidelberg
% Project:      TV restoration of spatially variant blur
%%
function A = forward_Sv1d(N, L, sigma, boundary)
    if ~is_positive_int(N)
        error("N must be a positive integer")
    end
    if ~is_positive_int(L)
        error("L must be a positive integer")
    end   
    if mod(L, 2) == 0
        error("L must be uneven")
    end
    if L > N
        error("N must be larger or equal L")
    end
    A = zeros(N, N+L-1);
    center = floor(L/2);

    % spatial location
    for si = 1:N
        % neighboring elements
        for ti = si-center:si+center
            % offset: ti + center = si + (ti - (si - center))
            A(si, ti+center) = kernel_Sv1d(si, ti, sigma);
        end
    end

    switch boundary
        case 'extended'
            % return full matrix

        case 'zero'
            A = A(:, center+1:N+L-1-center);

        case 'periodic'
            % split boundary columns from matrix A
            A_lb = A(:, 1:center);
            A_rb = A(:, end-center+1:end);
            A = A(:, center+1:N+L-1-center);
            
            % merge right boundary columns to first columns of A
            for col = 1:center
                idx_rb = A_rb(:, col) > 0;
                A(idx_rb, col) = A_rb(idx_rb, col);
            end

            % merge left boundary columns to last columns of A
            for col = N-center+1:N
                idx_lb = A_lb(:, col-(N-center)) > 0;
                A(idx_lb, col) = A_lb(idx_lb, col-(N-center));
            end

        case 'reflexive'
            A_lb_flip = flip(A(:, 1:center), 2);
            A_rb_flip = flip(A(:, end-center+1:end), 2);
            A = A(:, center+1:N+L-1-center);

            Z = zeros(N, size(A,2)-size(A_lb_flip,2)-size(A_rb_flip,2));
            A = A + horzcat(A_lb_flip, Z, A_rb_flip);

        otherwise
            error('Error: unsupported boundary condition')
    end
end

%%
function k = kernel_Sv1d(si, ti, sigma)
    c = 1 / (sqrt(2*pi) * sigma(si));
    k = c * exp(-0.5 * (ti-si)^2 / sigma(si)^2);
end

%%
function bool = is_positive_int(x)
    bool = isscalar(x) && x == floor(x) && x > 0;
end
