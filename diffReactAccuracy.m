%% define analytical solution
anal = @(x) ( (exp(3)/(exp(6)-1)) * (exp(3*x)-exp(-3*x)) );

%% get accuracy of estimated results

% set parameters
D      =  1;
lambda = -9;
f      =  0;

% test between 1 and 100 mesh elements (generates 10 numbers)
elementNumbers = round(logspace(0, 2, 10));
inaccuracies   = zeros(1,10); % empty array for inaccuracy scores

for i = 1:length(elementNumbers)
    
    % create mesh object
    mesh = OneDimLinearMeshGen(0,1,elementNumbers(i));
    
    % get estimated results
    results = finiteElementSolver(mesh,D,lambda,f,[1 0;mesh.ngn 1],[]);
    
    % generate 1000 x-points to test at
    xPoints = linspace(0,1,1000);
    
    % query analytical solution and results vector
    correct = anal(xPoints);
    attempt = interp1(mesh.nvec,results,xPoints);
    
    % get root of error squared
    error = sqrt((attempt - correct).^2);
    
    % save inaccuracy score
    inaccuracies(i) = sum(error) / 1000;
end

%% plot results
plot(elementNumbers, inaccuracies, '-x', 'Color', '#2e83dd');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('number of elements in mesh', 'FontSize', 12)
ylabel('inaccuracy as defined in report', 'FontSize', 12)
title('Inaccuracy of simulated result vs number of elements in mesh');