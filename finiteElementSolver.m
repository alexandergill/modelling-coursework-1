function cValues = finiteElementSolver(mesh, D, lambda, f, dirichletBCs, neumannBCs)
%finiteElementSolver Summary of this function goes here
%   Detailed explanation goes here

    %% initialise global matrix and source vector
    globalMatrix     = zeros(mesh.ngn);
    globalSrcVector  = zeros(mesh.ngn, 1); % source vector
    
    %% loop over mesh elements
    for i = 1:mesh.ne
        
        % calculate local Laplace matrix
        laplaceMatrix = LaplaceElemMatrix(D, i, mesh);
        % add it to global matrix
        globalMatrix(i:i+1, i:i+1) = globalMatrix(i:i+1, i:i+1)...
                                     + laplaceMatrix;
        
        % calculate local reaction matrix
        reactionMatrix = reactionElemMatrix(lambda, i, mesh);
        % subtract it from global matrix
        globalMatrix(i:i+1, i:i+1) = globalMatrix(i:i+1, i:i+1)...
                                     - reactionMatrix;
        
        % calculate local source vector and add to global vector
        sourceVector = sourceElemVector(f, i, mesh);
        globalSrcVector(i:i+1) = globalSrcVector(i:i+1) + sourceVector;
    end
    
    %% apply neumann boundary conditions
    % each condition is in the format [nID, grad], where grad is the known
    % gradient at the node nID
    for condition = neumannBCs.' % loop over rows
        
        % get required gradient and node ID
        nID  = condition(1);
        grad = condition(2);
        
        switch nID
            % BC is at the start of the mesh
            case 1
                globalSrcVector(nID) = globalSrcVector(nID) - grad;
            
            % BC is at the end of the mesh
            case mesh.ngn
                globalSrcVector(nID) = globalSrcVector(nID) + grad;
                
            otherwise
                error("can only enforce Neumann BC at mesh ends");
        end
    end
    
    %% apply dirichlet boundary conditions
    % each condition is in the format [nID, c], where c is the known
    % concentration at the node nID
    for condition = dirichletBCs.'
        
        % get required concentration and node ID
        nID = condition(1);
        c   = condition(2);
        
        % get the identity matrix for the number of nodes
        I = eye(mesh.ne + 1);
        
        % set corresponding matrix row to a slice of the identity matrix
        globalMatrix(nID, :) = I(nID, :);
        
        % over-write the corresponding element in the globalVector
        globalSrcVector(nID) = c;
    end
    
    %% calculate c vector using c = M^-1 f
    cValues = globalMatrix \ globalSrcVector;
end
