function reactionMatrix = reactionElemMatrix(lambda, eID, msh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % verify that element eID is in the mesh
    if eID > msh.ne
        error('mesh element ID is not in mesh')
    end
    
    % get Jacobian for this element
    J = msh.elem(eID).J;
    
    % create empty matrix
    reactionMatrix = zeros(2);
    
    % fill the matrix
    for m = 0:1
        for n = 0:1
            element = lambda * get_dPsi_by_dXi(n) * ...
                      get_dPsi_by_dXi(m) * J;
            
            % put element in matrix
            i = m + 1; j = n + 1;
            reactionMatrix(i, j) = element;
        end
    end
end
