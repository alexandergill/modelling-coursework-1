function sourceVector = sourceElemVector(f, eID, msh)

    % verify that element eID is in the mesh
    if eID > msh.ne
        error('mesh element ID is not in mesh')
    end
    
    % get Jacobian for this element
    J = msh.elem(eID).J;
    
    % create 2x1 source vector with both elements = fJ
    sourceVector = ones(2, 1) * f * J;
    