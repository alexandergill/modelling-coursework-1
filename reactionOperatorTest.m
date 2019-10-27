%% Test 1: test symmetry of the matrix
% Test that this matrix is symmetric
tol = 1e-14;
eID = 1; % element ID
lambda = 10; % reaction coefficient
msh = OneDimLinearMeshGen(0,1,10);

reactionMatrix = reactionElemMatrix(lambda, eID, msh);

assert(abs(reactionMatrix(1,2) - reactionMatrix(2,1)) <= tol)
