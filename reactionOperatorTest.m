%% Test 1: test symmetry of the matrix
% Test that this matrix is symmetric
tol = 1e-14;
eID = 1; % element ID
lambda = 10; % reaction coefficient
msh = OneDimLinearMeshGen(0,1,10);

reactionMatrix = reactionElemMatrix(lambda, eID, msh);

assert(abs(reactionMatrix(1,2) - reactionMatrix(2,1)) <= tol)

%% Test 2: test that one matrix is evaluted correctly
% Test that a known matrix is calculated correctly
tol = 1e-14;
eID = 1; % element ID
lambda = 10; % reaction coefficient
msh = OneDimLinearMeshGen(0,1,3);

attempt = reactionElemMatrix(lambda, eID, msh);

correct = [lambda / 9 , lambda / 18; ...
           lambda / 18, lambda / 9];

diff = attempt - correct;
diffnorm = sum(sum(diff.*diff)); % total squared error
assert(abs(diffnorm) <= tol)