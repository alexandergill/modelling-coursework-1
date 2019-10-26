function [dPsi_by_dXi] = get_dPsi_by_dXi(n)
%GET_DPSI_BY_DXI Get the differential of basis function Psi_n

    switch n
        case 0
            dPsi_by_dXi = -0.5;
        case 1
            dPsi_by_dXi = 0.5;
        otherwise
            error('invalid basis function')
    end
end

