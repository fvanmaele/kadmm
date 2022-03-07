%% forward_Sv2d
% Forward operator for 2d separable Gaussian blur, discretized to size LxL
%
%% Input
% M:            amount of rows in image
% N:            amount of columns in image
% L:            kernel size
% gamma:        blur parameter
% boundary:     type of boundary condition
% blurType:     'corners', 'center', 'gradient' or 'constant'
%
%% Output
% A:            forward operator
%
%% References
% File:         forward_Sv2d.m
% Author:       Ferdinand Vanmaele, University of Heidelberg
% Project:      TV restoration of spatially variant blur
%%
function A = forward_Sv2d(M, N, L, gamma, boundary, blurType)
    % variance depending on spatial location
    switch blurType
        case 'corners'
            sigma1 = @(s) gamma * abs(0.5 - s/M) + 0.5;
            sigma2 = @(s) gamma * abs(0.5 - s/N) + 0.5;
            
        case 'center'
            sigma1 = @(s) -gamma * abs(0.5 - s/M) + 2.5;
            sigma2 = @(s) -gamma * abs(0.5 - s/N) + 2.5;
            
        case 'gradient'
            sigma1 = @(s) gamma * 0.5 * (0.5 - s/M) + 1.25;
            sigma2 = @(s) gamma * 0.5 * (0.5 - s/N) + 1.25;
            
        case 'constant'
            sigma1 = @(s) gamma;
            sigma2 = @(s) gamma;

        otherwise
            error("blur type not supported")
    end

    A1 = sparse(forward_Sv1d(M, L, sigma1, boundary));
    A2 = sparse(forward_Sv1d(N, L, sigma2, boundary));
    
    % kronecker product for seperable blur
    A = kron(A2, A1);
end