function reactionMatrix = reactionElemMatrix(lambda, eID, msh)
%UNTITLED2 Generates reaction element matrix
%   Generates the reaction element matrix for element 'eID' of mesh 'msh'
%   given a reaction coefficient 'lambda'

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
            
            % get basis functions
            psi_m = getPsi(m);
            psi_n = getPsi(n);
            
            % construct function and calculate integral
            func = @(xi) lambda .* psi_m(xi) .* psi_n(xi) .* J;
            elem = integral(func, -1, 1);
            
            % put integral in matrix
            i = m + 1; j = n + 1;
            reactionMatrix(i, j) = elem;
        end
    end
end
