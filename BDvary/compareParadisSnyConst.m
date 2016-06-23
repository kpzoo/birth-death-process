% Function to simulate and compare the Snyder constant rate case to the
% Paradis fitting scheme to assess accuracy, uses Tanja method
function [est, estParadis, qmarg] = compareParadisSnyConst()

% Assumptions and modifications
% - includes comparison of Paradis cdf estimates for constant rate
% - condition on a fixed time, T and n (previously drawed n)
% - sampling prob of 1
% - discretised time and used T, mu, a, b (= k in paper) from Nee 1994


%% Initial parameters and definitions for constant rate birth-death

% Parameter names and true values (assume lam > mu)
lam = 0.1;
mu = 0.5*lam;
x = [lam mu];
xdefs = {'lam', 'mu'};
sig = lam - mu;
rho = mu/lam;

% Points for discretised inference space
numRV = length(x);
mi = 30*ones(size(x));

% Space symmetrically centred on true values
xmax = 10*x;
xmin = 0.1*x;

% Calculate parameter space
m = prod(mi);
xset = cell(1, 1);
for i = 1:numRV
    xset{i} = linspace(xmin(i), xmax(i), mi(i));
end

% Define no. tips
n = 200;
nLin = 1:n;
nData = n-1;


%% Simulate a true birth-death reconstructed tree from Hartmann 2008

% Simulate the constant BD process according to Hartmann 2008
[T, tbranch] = genBirDeaConst(n, sig, rho);

% Convert branch times into coalescents (0 is present)
if max(tbranch) >= T
    error('The branching times are not proper');
end
tcoal = sort([0 tbranch]);

% Convert times forward from T
tspec = sort(T - tcoal);
tspec = tspec - tspec(1); % set zero as first speciation
Tsp = max(tspec);

% Discretise time points for numerical integration
nVals = 20000;
t = linspace(0, Tsp, nVals);

% Calculate Paradis 2011 theoretical cdf
Dfn = @(lam, mu, t, T) exp((lam - mu)*t).*(lam - mu*exp(-(lam - mu)*(T*ones(size(t)) - t))).^(-2);
tDden = trapz(t, Dfn(lam, mu, t, Tsp), 2);
tDnum = cumtrapz(t, Dfn(lam, mu, t, Tsp), 2);
tDcdf = tDnum/tDden;
    
% Calculate a Paradis estimate by brute force over same prior range
nPts = 300;
estParadis = testParadis(xmin, xmax, nPts, n, tspec, Dfn, t, Tsp, xdefs, x, tDcdf);
disp('Paradis estimate');
disp(estParadis);



%% Snyder filtering of the reconstructed process

% Set uniform joint prior q0 and maximum speciation time
q0 = ones(1, m)/m;

% Create a matrix of identifiers to tell which xset{i} values are used for
% each entry in N(t) and lam(t) calculations
IDMx = zeros(numRV, m);
% Initialise with first variable which has no element repetitions
idxset = 1:mi(1);
IDMx(1, :) = repmat(idxset, 1, m/mi(1));
for i = 2:numRV
    % For further variables numReps gives the number of set repetitions
    % while kronVec gives the number of element repetitions
    idxset = 1:mi(i);
    numReps = m/prod(mi(1:i));
    kronVec = ones(1, prod(mi(1:i-1)));
    IDMx(i, :) = repmat(kron(idxset, kronVec), 1, numReps);
end

% Get the values corresponding to the matrix
xsetMx = zeros(numRV, m);
for i = 1:numRV
    xsetMx(i, :) = xset{i}(IDMx(i, :));
end

% Use known analytic constant rate solution
[est, qmarg] = testConst(x, tspec, m, nLin, nData, xset, q0,...
    xsetMx, numRV, Tsp, xdefs, IDMx, mi);
    