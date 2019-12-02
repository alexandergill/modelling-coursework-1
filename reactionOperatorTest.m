%% Test 1: test symmetry of the matrix
% Test that this matrix is symmetric
tolerance = 1e-14;
eID       = 1; % element ID
lambda    = 10; % reaction coefficient
mesh      = OneDimLinearMeshGen(0,1,10);

reactionMatrix = reactionElemMatrix(lambda, eID, mesh);

assert(abs(reactionMatrix(1,2) - reactionMatrix(2,1)) <= tolerance)

%% Test 2: test that one matrix is evaluted correctly
% Test that a known matrix is calculated correctly
tolerance = 1e-14;
eID       = 1; % element ID
lambda    = 10; % reaction coefficient
mesh      = OneDimLinearMeshGen(0,1,3);

attempt = reactionElemMatrix(lambda, eID, mesh);

correct = [lambda / 9 , lambda / 18; ...
           lambda / 18, lambda / 9];

diff = attempt - correct;
diffnorm = sum(sum(diff.*diff)); % total squared error
assert(abs(diffnorm) <= tolerance)

%% Test 3: test that the same matrix is produced for same-sized elements
tolerance = 1e-14;
eID       = 2;
lambda    = 9;
mesh      = OneDimLinearMeshGen(0,1,4);

matrix2 = reactionElemMatrix(lambda, eID, mesh);

eID = 3;
matrix3 = reactionElemMatrix(lambda, eID, mesh);

% check each element of the matrix is the same
for i=1:2
    for j=1:2
        assert(abs(matrix2(i,j)-matrix3(i,j)) <= tolerance)
    end
end