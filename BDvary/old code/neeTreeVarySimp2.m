% Script to draw some speciation times from using the time varying rate
% algorithms from Hohna 2013 but simplified to a constant rate process

% Assumptions and modifications
% - modified to allow a simpler lam(t) function
% - extended neeTreeConstSimp2 with Nee for time varying 3 param case
% - sampling prob of 1 and constant death rate
% - discretised time and used T, mu, a, b (= k in paper) from Nee 1994

clearvars
clc
close all

tic;

%% Basic parameters and functions

% Boolean to control functions
func = 1;
funcStr = {'Nee', 'linear'};
funlab = funcStr(func);

% Pick parameters based on functions and define the functions
switch(func)
    case 1
        % Nee 1994 variable birth, constant death
        mi = [10 10 10];
        minSpace = [0.01 0.01 0.01];
        maxSpace = [10 50 0.1];
        
        % Set true parameters from Nee 1994
        mu = 0.075;
        b = 4;
        a = 0.2;
        x = [a b mu];
        xdefs = {'a', 'b', 'mu'};
        
        % Define net mu-lam rate function from Nee 1994 and its integral
        lamb = @(x, tx) x(2)*x(3)./(1 + x(1)*tx);
        mub = @(x, tx) x(3)*ones(size(tx));
        netB = @(x, tx) mub(x, tx) - lamb(x, tx);
        intB = @(x, t1, t2) x(3)*(t2 - t1) + (x(3)*x(2)/x(1))*(log((x(1)*t1 + 1)./(x(1)*t2 + 1)));
        
        
    case 2
        % Linear birth, constant death
        mi = [30 30];
        minSpace = [0.01 0.1];
        maxSpace = [0.1 1];
        xdefs = {'a', 'mu'};
        
        % Define rate functions and integrals
        lamb = @(x, tx) x(1)*tx + x(2);
        mub = @(x, tx) x(2)*ones(size(tx));
        netB = @(x, tx) mub(x, tx) - lamb(x, tx);
        intB = @(x, t1, t2) - 0.5*x(1)*(t2.^2 - t1.^2);
        
    otherwise
        error('Incorrect function type parsed');
end

% Initialisation of sim parameters
numRV = length(mi);
m = prod(mi);
n = 100;

% Generate space and true values of parameters if needed
xset = cell(1, 1);
if ~exist('x', 'var')
    x = zeros(1, numRV);
    medianX = 1;
else
    medianX = 0;
end
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
    if medianX
        % Assign true value if not done already as median of range
        x(i) = median(xset{i});
    end
end

% Constant rate function r(t, s) = e^(int mu-lam) and integrand of PtT
rst = @(t1, t2, x) exp(intB(x, t1, t2));
intT = @(t1, t2, x) lamb(x, t1).*rst(t1, t2, x);

%% Simulate a true birth-death tree

% Generate a suitable T assuming uniform prior across Tset
nT = 5000;
Tset = linspace(0, 100, nT);
PTn = zeros(size(Tset));

for i = 1:nT
    % For each possible T get the PTn entry
    Ti = Tset(i);
    P0Ts = (1 + integral(@(t2)intT(0*ones(size(t2)), t2, x), 0, Ti))^(-1);
    % Get the P(T|n1) entry for each T (unnormalised)
    PTn(i) = ((1 - P0Ts*rst(0, Ti, x)).^(n-1)).*(P0Ts^2)*rst(0, Ti, x);
end

% Take T as maximum value (don't need to normalise)
T = Tset(PTn == max(PTn));

% Discretise time points for numerical integration
nVals = 20000;
t = linspace(0, T, nVals);
PtT = zeros(size(t));
ertT = zeros(size(t));
for i = 1:nVals
    PtT(i) = (1 + integral(@(t2)intT(t(i), t2, x), t(i), T))^(-1);
    ertT(i) = rst(t(i), T, x);
end

% Get Pr(N(T) = 1 | N(t) = 1) = p1tT
p1tT = (PtT.^2).*ertT;

% Probability of n lineages at T given 1 at t = 0
nLinSet = 1:500;
pn0T = zeros(size(nLinSet));
P0T = PtT(1);
er0T = ertT(1);
for i = nLinSet
    pn0T(i) = ((1 - P0T*er0T).^(i-1))*(P0T.^2)*er0T;
end

%% Hohna 2013 method for drawing speciation times

% Using discretised form of eq 9 from Hohna 2013 where t^* replaced by
% interpolated t that satisfies 9 - take n-1 draws
r = rand(1, n-1);
tden = trapz(t, lamb(x, t).*p1tT, 2);
tnum = cumtrapz(t, lamb(x, t).*p1tT, 2);
tcdf = tnum./tden;

% Use interpolation to find speciation times corresponding to rand values
tspec = interp1(tcdf, t, r);
tspec = sort([0 tspec]);

% Remove duplicate times and set the data length
tspec = unique(tspec);
n1 = length(tspec);
if n ~= n1
    disp(['Change in n to ' num2str(n1) ' from ' num2str(n)]);
    n = n1;
end
nLin = 1:n;
nData = n-1;

% tspec = tspec(2:end);
% n = n-1;
% nLin = 2:n;
% nData = n-1;


%% Perform Snyder estimation for 3 parameters using Nee 1994 Poisson rate

% Set uniform joint prior q0 and maximum speciation time
q0 = ones(1, m)/m;
Tsp = max(tspec);

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
    [tsol, qsol] = ode113(@(ts, y) odeSnyIntVaryFn(ts, y, xsetMx, numRV, Tsp, nLinCurr, intT, lamb),...
        [tspec(i) tspec(i+1)], qev(i, :)', options);
    
    % Assign the output values of time and posteriors
    qset{i} = qsol;
    tset{i} = tsol;
    elemLen(i) = length(tsol);
    
    % Perturb the q posterior for the new event
    lampert = getNeeTimeVaryIntFn(tsol(end), xsetMx, numRV, Tsp, nLinCurr, intT, lamb);
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :).*lampert./(qev(i+1, :)*lampert');
    
    disp(['Finished: ' num2str(i) ' of ' num2str(nData)]);
    
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

% Display structure for parameters
est.x = x;
est.xhat = xhat(end, :);
est.perc = 100*(1 - est.xhat./x);
est.lab = xdefs;
disp(est);

% Final posterior and marginalisations
qnlast = qn(end, :);
[qmarg, ~] = marginalise(numRV, IDMx, qnlast, mi);

% True birth and death rates
mut = mub(x, tn);
lamt = lamb(x, tn);

% Estimated rates based on function type
switch(func)
    case 1
        % Nee 1994 rates
        muhat = est.xhat(3)*ones(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) = qnlast*(xsetMx(2, :).*xsetMx(3, :)./(1 + xsetMx(1, :)*tn(i)))';
        end
    case 2
        % Linear birth rate, constant death
        muhat = est.xhat(2)*ones(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) = qnlast*(xsetMx(1, :)*tn(i) + xsetMx(2, :))';
        end
end

% Collect rates for plotting
rates = [lamt lamhat mut muhat];
maxr = max(max(rates));
minr = min(min(rates));

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


% % Plot the true and estimated rate differences
% figure;
% subplot(2, 1, 1);
% plot(tn, lamt - mut, tn, lamhat - muhat);
% ylabel('\lambda - \mu');
% legend('real', 'estimated', 'location', 'best');
% xlabel('time');
% title('Rate differences');
% grid;
% subplot(2, 1, 2);
% plot(tn, cumtrapz(tn, lamt - mut, 1), tn, cumtrapz(tn, lamhat - muhat, 1));
% ylabel('integral of \lambda - \mu');
% legend('real', 'estimated', 'location', 'best');
% xlabel('time');
% title('Rate integral differences');
% grid;

% Plot characteristics of the no. lineages simulated at T
figure;
subplot(2, 1, 1);
plot(nLinSet, pn0T);
xlabel('no. lineages at T, n');
ylabel('Pr(N(T) = n | N(0) = 1)');
title(['Probability on no. lineages at T = ' num2str(T)]);
grid;
subplot(2, 1, 2);
stairs(tspec, 1:n);
xlabel('time');
ylabel('no lineages');
title('Lineages through time plot');
grid;

% Plot birth and death rates
figure;
[hAx, hLine1, hLine2] = plotyy(tn, [lamt lamhat], tn, [mut muhat]);
xlabel('time');
xlim([tn(1) tn(end)]);
ylabel(hAx(1),'birth rate') % left y-axis
ylabel(hAx(2),'death rate') % right y-axis

% Birth rates
hLine1(1).LineStyle = '-';
hLine1(1).Color = 'b';
hLine1(2).LineStyle = '-';
hLine1(2).Color = 'r';

% Death rates
hLine2(1).LineStyle = '-.';
hLine2(1).Color = 'b';
hLine2(2).LineStyle = '-.';
hLine2(2).Color = 'r';

legend('\lambda', '\lambda est', '\mu', '\mu est', 'location', 'best');
title('Estimated rates');
grid;