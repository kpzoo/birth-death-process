% Script to draw some speciation times from a rate time varying birth-death
% process according to the example in Nee 1994 where deaths are constant

% Assumptions and modifications
% - NOT WORKING
% - modified to use cumtrapz
% - sampling prob of 1 and constant death rate
% - discretised time and used T, mu, a, b (= k in paper) from Nee 1994

clearvars
clc
close all

tic;

%% Simulate a true tree

% Set true parameters from Nee 1994
T = 300;
mu = 0.075;
b = 4;
a = 0.2;
n1 = 100;

% Discretise time points for numerical integration
nVals = 20000;
t = linspace(0, T, nVals);

% Time varying birth rate, lam
lam = b*mu./(1 + a*t);

% P(t, T) with 1st entry as P(0, T)
PtT = getPtT(lam, mu*ones(size(t)), t);
%PtTd = getPtTDense(a, b, mu, t);

% Calculate r0t then rtT fom r0T = r0t + rtT
%r0t = cumtrapz(t, mu - lam, 2);

% Iwasa
r0t = cumtrapz(t, lam - mu, 2);

rtT = r0t(end) - r0t;

% Get Pr(N(T) = 1 | N(t) = 1) = p1tT
p1tT = PtT.*PtT.*exp(rtT);

% Avg no. lineages from Nee 1994
navg = exp(-r0t).*PtT/PtT(1);
Mavg = 1 + cumtrapz(t, exp(-r0t).*lam);
nLinT = [navg(end) Mavg(end)];

% Probability of n lineages at T given 1 at t = 0
nLinSet = 1:100;
pn0T = zeros(size(nLinSet));
P0T = PtT(1);
r0T = rtT(1);
for i = nLinSet
    pn0T(i) = ((1 - P0T*exp(r0T)).^(i-1))*(P0T.^2)*exp(r0T);
end

% Using discretised form of eq 9 from Hohna 2013 where t^* replaced by
% nearest element of t that satisfies 9 - take n-1 draws
r = rand(1, n1-1);
% Cumulative integral then find t for each integ: cumtrapz <= integ
svals = cumtrapz(t, lam.*p1tT, 2);
svals = svals/svals(end);
sset = zeros(1, n1-1);
for i = 1:n1-1
    sset(i) = t(find(svals >= r(i), 1, 'first'));
end

% Speciation times for n species
tspec = sort([0 sset]);

% Remove duplicate times
tspec = unique(tspec);
n = length(tspec);
nLin = 1:n;
nData = n-1;

% Obtain lineages and times across t vector
nt = zeros(size(t));
for i = 1:n
    nt(tspec(i) == t) = 1;
end
nt = cumsum(nt);


% Calculate mean no. lineages at any t from inhomogeneous Poisson rate
pRate = lam.*PtT.*nt;
Lamt = cumtrapz(t, pRate);

%% Perform Snyder estimation for 3 parameters using Nee 1994 Poisson rate

% Initialisation of sim parameters
mi = [10 10 10];
numRV = length(mi);
m = prod(mi);

% Parameter search space
minSpace = [0.01 0.1 0.001];
maxSpace = [1 10 0.1];

% minSpace(1) = a;
% maxSpace(1) = a;
% minSpace(2) = b;
% maxSpace(2) = b;
% mi = [1 1 100];
% numRV = length(mi);
% m = prod(mi);


xdefs = {'a', 'b', 'mu'};
x = [a b mu];
xset = cell(1, 1);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
end

% Set uniform joint prior q0 and maximum speciation time
q0 = ones(1, m)/m;
Tsp = max(tspec);
%Tsp = T;

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
    [tsol, qsol] = ode113(@(ts, y) odeSnyBDvary(ts, y, xsetMx, numRV, Tsp, nLinCurr),...
        [tspec(i) tspec(i+1)], qev(i, :)', options);
    
    % Assign the output values of time and posteriors
    qset{i} = qsol;
    tset{i} = tsol;
    elemLen(i) = length(tsol);
    
    % Perturb the q posterior for the new event
    lampert = getNeeTimeVaryRate(tsol(end), xsetMx, numRV, Tsp, nLinCurr);
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :).*lampert./(qev(i+1, :)*lampert');
    
    disp(num2str(i));
    
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

% Simulation time
tsim = toc/60;
disp(['Simulation time = ' num2str(tsim) ' mins']);

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

% Estimated and true birth and death rates
mut = mu*ones(size(tn));
lamt = x(2)*x(3)./(1 + x(1)*tn);
muhat = est.xhat(3)*ones(size(tn));
lamhat = zeros(size(lamt));
for i = 1:length(tn)
    fac = xsetMx(2, :).*xsetMx(3, :)./(1 + tn(i)*xsetMx(1, :));
    lamhat(i) = qn(i, :)*fac';
end
lamhat = reshape(lamhat, size(lamt)); % make sure it is column
rates = [lamt lamhat mut muhat];
maxr = max(max(rates));
minr = min(min(rates));

% Plot of birth-death rates
figure;
plot(tn, lamt, tn, lamhat, tn, mut, tn, muhat);
ylim([minr maxr]);
xlim([min(tn) max(tn)]);
grid;
xlabel('time');
ylabel('birth and death rates');
legend('\lambda (t)', '\lambda (t) est', '\mu (t)', '\mu(t) est', 'location', 'best');

% Plot the marginal posteriors
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

% Plot the true and estimated rate differences
figure;
subplot(2, 1, 1);
plot(tn, lamt - mut, tn, lamhat - muhat);
ylabel('\lambda - \mu');
legend('real', 'estimated', 'location', 'best');
xlabel('time');
title('Rate differences');
grid;
subplot(2, 1, 2);
plot(tn, cumtrapz(tn, lamt - mut, 1), tn, cumtrapz(tn, lamhat - muhat, 1));
ylabel('integra of \lambda - \mu');
legend('real', 'estimated', 'location', 'best');
xlabel('time');
title('Rate integral differences');
grid;

% Plot characteristics of the no. lineages simulated at T
figure;
subplot(2, 1, 1);
plot(nLinSet, pn0T);
xlabel('no. lineages at T, n');
ylabel('Pr(N(T) = n | N(0) = 1)');
title(['Probability on no. lineages at T = ' num2str(T)]);
grid;
subplot(2, 1, 2);
plot(t, navg, t, Mavg, t, Lamt, t, nt);
xlabel('time');
ylabel('avg. no lineages');
legend('Nee avg', 'Kendall avg', 'inhomo avg', 'sim no.', 'location', 'best');
title('Different averages on lineages across time');
grid;