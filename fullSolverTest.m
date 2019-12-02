%% Test 1: test solution to laplace's equation with dirichlet BCs
tolerance = 1e-14; % tolerance for comparing floating-point values
D         = 1;     % diffusion coefficient
mesh      = OneDimLinearMeshGen(0,1,4); % create mesh object

% get results
c = finiteElementSolver(mesh,1,0,0,[1 2; mesh.ngn 0],[]);

% check against analytical solution
anal = @(x) 2*(1-x); % analytical solution
for i = 1:mesh.ngn   % loop over mesh nodes
    
    % check solutions are the same
    assert(abs(c(i)-anal(mesh.nvec(i)))<=tolerance); 
end

%% Test 2: test neumann boundary condition at x=0
tolerance       = 1e-14; % tolerance for comparing floating-point values
D               = 1;     % diffusion coefficient
mesh            = OneDimLinearMeshGen(0,1,4); % create mesh object
correctGradient = 2;     % required neumann BC

% get results
c = finiteElementSolver(mesh,1,0,0,[mesh.ngn 0],[1 correctGradient]);

% get gradient at x=0
gradient = (c(2) - c(1)) / (mesh.nvec(2) - mesh.nvec(1));
assert(abs(gradient - correctGradient) <= tolerance);

%% Test 3: test neumann boundary condition at x=1
tolerance       = 1e-14; % tolerance for comparing floating-point values
D               = 1;     % diffusion coefficient
mesh            = OneDimLinearMeshGen(0,1,4); % create mesh object
correctGradient = 2;     % required neumann BC

% get results
c = finiteElementSolver(mesh,1,0,0,[1 0],[mesh.ngn correctGradient]);

% get gradient at x=0
gradient = (c(2) - c(1)) / (mesh.nvec(2) - mesh.nvec(1));
assert(abs(gradient - correctGradient) <= tolerance);

%% Test 4: test source vector contains two values
tolerance = 1e-14; % tolerance for comparing floating-point values
f         = 9;     % source coefficient
mesh      = OneDimLinearMeshGen(0,1,4); % create mesh object

sourceVector = sourceElemVector(f,1,mesh);

assert(length(sourceVector) == 2);

%% Test 5: test source vector elements are equal
tolerance = 1e-14; % tolerance for comparing floating-point values
f         = 9;     % source coefficient
mesh      = OneDimLinearMeshGen(0,1,4); % create mesh object

sourceVector = sourceElemVector(f,1,mesh);

assert(abs(sourceVector(1)-sourceVector(2)) <= tolerance);

%% Test 6: test that source vector values are correct for known solution
tolerance = 1e-14; % tolerance for comparing floating-point values
f         = 9;     % source coefficient
mesh      = OneDimLinearMeshGen(0,1,4); % create mesh object

sourceVector = sourceElemVector(f,1,mesh);
knownVector  = [1.125; 1.125];

for i = 1:length(sourceVector)
    assert(abs(sourceVector(i) - knownVector(i)) <= tolerance);
end

%% Test 7: test that d_psi_by_d_xi gives the correct values
tolerance = 1e-14; % tolerance for comparing floating-point values

% get values
d_psi_by_d_xi_0 = get_dPsi_by_dXi(0);
d_psi_by_d_xi_1 = get_dPsi_by_dXi(1);

% check they're correct
assert(abs(d_psi_by_d_xi_0 + 0.5) <= tolerance);
assert(abs(d_psi_by_d_xi_1 - 0.5) <= tolerance);