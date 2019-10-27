function [psi] = getPsi(n)
%GETPSI Summary of this function goes here
%   Detailed explanation goes here
    switch n
        case 0
            psi = @(xi) (1 - xi) / 2;
        case 1
            psi = @(xi) (1 + xi) / 2;
        otherwise
            error('psi%d is not a valid function', n)
end

