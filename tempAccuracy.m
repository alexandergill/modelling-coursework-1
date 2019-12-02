%% get estimated gradients

% set parameters
k     = 1.0100e-05;
Q     = 1.5;
T_L   = 294.15;
T_in  = 293.15;
T_out = 323.15;

% test between 1 and 1000 mesh elements (generates 10 numbers)
elementNumbers = round(logspace(0, 3, 10));
results        = zeros(1,10); % empty array for gradient values

for i = 1:length(elementNumbers)
    
    % create mesh object
    mesh = OneDimLinearMeshGen(0,0.01,elementNumbers(i));
    
    % get estimated gradient at x = 0 (node 1)
    tempDist = finiteElementSolver(mesh, k, -Q, Q * T_L, ...
                                   [1 T_out; mesh.ngn T_in],[]);
    gradient = (tempDist(2) - tempDist(1)) / (mesh.nvec(2) - mesh.nvec(1));
    results(i) = gradient;
end

%% plot results
% /1000 to get K/mm
plot(elementNumbers, results/1000, '-x', 'Color', '#2e83dd');
set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
xlabel('number of elements in mesh', 'FontSize', 12)
ylabel('gradient at x=0 / Kmm^{-1}', 'FontSize', 12)
title('Initial gradient vs number of elements in mesh', 'FontSize', 16);