function [diffusionMatrix] = LaplaceElemMatrix(D,eID,msh)
%LaplaceMatrix Generates local element matrix for the diffusion operator
%   Returns a 2x2 matrix TODO:

    % verify that element eID is in the mesh
    if eID > msh.ne
        error('mesh element ID is not in mesh')
    end
    
    % get Jacobian for this element
    J = msh.elem(eID).J;
    dXi_by_dx = 1 / J;
    
    % create empty matrix
    diffusionMatrix = zeros(2);
    
    % fill empty matrix
    
    % m is the row index
    for m = 0:1
        % n is the column index
        for n = 0:1
            % calculate element
            element = D * get_dPsi_by_dXi(n) * ...
                      dXi_by_dx * get_dPsi_by_dXi(m) * ...
                      dXi_by_dx * J * 2;
            %                         ^
            % x2 to perform integral from -1 to 1
            
            % hack because ShatLab doesn't allow 0 indices
            i = m+1;
            j = n+1;
            
            % shove it in the matrix
            diffusionMatrix(i,j) = element;
        end
    end
end

