function [psi] = getPsi(n)
%GETPSI Returns the basis function psi_n where n = {0, 1}
    switch n
        case 0
            psi = @(xi) (1 - xi) / 2;
        case 1
            psi = @(xi) (1 + xi) / 2;
        otherwise
            error('psi%d is not a valid function', n)
end

