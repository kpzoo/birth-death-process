% Script to simulate the birth-death reconstructed process from the time
% varying birth perocess of Nee et al 1994

% Assumptions and modifications
% - constant birth and death rates and sampling of 1


clc
close all
clearvars

%% Initialisation and simulation of birth-death process

% Initialise Snyder filtering problem for estimation of [rho sig]
n = 100;
mi = [100 100];
numRV = length(mi);
m = prod(mi);
nData = n-1;
minSpace = [0.01 0.01];
maxSpace = [0.99 100];
xdefs = {'rho', 'sigma'};

% Get the instance of the random variables for rate function, lam
xset = cell(1, 1);
x = zeros(1, numRV);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
    %x(i) = datasample(xset{i}, 1);
    x(i) = xset{i}(round(mi(i)/2));
end

% Set parameters for true constant birth-death 
rho = x(1);
sig = x(2);
lam = sig/(1 - rho);
mu = lam - sig;

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
nLin = 1:n;
Tsp = max(tspec); 

% Check of this time
tcheck = sort(max(tcoal) - tcoal);
if any(abs(tcheck - tspec) > 10^-(10))
    disp('Times are inconsistent in forward time');
end

%% Snyder causal filtering

% Set uniform prior q0
q0 = ones(1, m)/m;

% Posterior vectors on events
qev = zeros(nData+1, m);
qev(1, :) = q0;
tev = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);
options = odeset('NonNegative', 1:m);
elemLen = zeros(1, nData);

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


% Loop across coalescent events and perform filtering 
for i = 1:nData
    % Current no. lineages
    nLinCurr = nLin(i);
    
    % Solve linear ODEs continuously with setting of options, no Q
    [tsol, qsol] = ode113(@(ts, y) odeSnyBDconst(ts, y, xsetMx, numRV, Tsp, nLinCurr),...
        [tspec(i) tspec(i+1)], qev(i, :)', options);
    
    % Assign the output values of time and posteriors
    qset{i} = qsol;
    tset{i} = tsol;
    elemLen(i) = length(tsol);
    
    % Perturb the q posterior for the new event
    lampert = getNeeRate(tsol(end), xsetMx, numRV, Tsp, nLinCurr);
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :).*lampert./(qev(i+1, :)*lampert');
    
end

% Get full length of ODE solution data and assign appending vectors
lenFull = sum(elemLen);
stop = 0;
qn = -ones(lenFull, m);
tn = -ones(lenFull, 1);

% Append the cell based posterior and time data into a single structure
for i = 1:nData
    % Loop calculates the start and end points along the array successively
    % and then assigns the appropriate cell element
    start = stop + 1;
    stop = elemLen(i) + start - 1;
    tn(start:stop) = tset{i};
    qn(start:stop, :) = qset{i};
end

% Time range
tnmin = min(tn);
tnmax = max(tn);
tnlen = length(tn);
    
%% Analyse and plot results

% Conditional mean
xhat = qn*xsetMx';
xhatmean = xhat;
% Conditional variance
xhatvar = qn*(xsetMx'.^2) - xhat.^2;
xhatstd = sqrt(xhatvar);
xstdub = xhatmean + 2*xhatstd;
xstdlb = xhatmean - 2*xhatstd;

% Display structure for parameters
est.x = x;
est.xhat = xhat(end, :);
est.lab = xdefs;

% Final posterior and marginalisations
qnlast = qn(end, :);
[qmarg, probSums] = marginalise(numRV, IDMx, qnlast, mi);

% Define space then multiply by posterior to estimate birth and death rates
birdea.label = {'lam', 'mu'};
birdea.real = [lam mu];
lamq = xsetMx(2, :)./(1 - xsetMx(1, :));
muq = lamq - xsetMx(2, :);
lamq = qnlast*lamq';
muq = qnlast*muq';
birdea.est = [lamq muq];
birdea.perc = 100*(1 - birdea.est./birdea.real).^2;
disp('Real and estimated rates are:');
disp(birdea);

% Plot joint bivariate distribution
if numRV == 2
    qtemp = reshape(qnlast, mi);
    if ~isnan(qtemp)
        figure;
        meshz(xset{1}, xset{2}, qtemp');
        xlabel('\rho');
        ylabel('\sigma');
        zlabel('P(\rho, \sigma | data)');
        title(['Joint posterior with true [\rho \sigma] = [' num2str(x(1)) ' ' num2str(x(2)) ']']);
    else
        warning('Mat:qNaN', 'qtemp has NaN values');
    end
end

% Plot marginals for each parameter and real value
figure;
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    plot(xset{i}, qmarg{i}, 'b', [x(i) x(i)], [0 max(qmarg{i})], 'k',...
        [xhat(end, i) xhat(end, i)], [0 max(qmarg{i})], 'r');
    xlabel(['x_' num2str(i)]);
    ylabel(['P(x_' num2str(i) ')']);
    legend('marginal', 'true', 'cond mean', 'location', 'best');
    grid;
end

% The parameters estimated and the std deviation bounds
figure;
for i = 1:numRV
    subplot(numRV, 1, i);
    plot(tn, x(i)*ones(size(tn)), 'r', tn, xhatmean(:, i), 'k', tn, xstdlb(:, i),...
        'b', tn, xstdub(:, i), 'b');
    xlabel('time');
    ylabel(['x_' num2str(i)]);
    title(['Estimate of x_' num2str(i) 'for [n m] = [' num2str(n) ' ' num2str(m) ']']);
    legend('true', 'estimate', '-2 std', '+2 std', 'location', 'best');
    xlim([tnmin tnmax]);
end
