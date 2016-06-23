% Script to simulate and estimate time varying birth-death models with
% comparison to analytic results in the constant rate case

% Assumptions and modifications
% - condition on a fixed time, T instead of n
% - uses algorithms from Hohna 2013
% - sampling prob of 1
% - discretised time and used T, mu, a, b (= k in paper) from Nee 1994

clearvars
clc
close all

tic;

% Booleans to control performance
nBool = 1; % just set n instead of drawing it

%% Initial parameters and definitions

% Define rate function type
rateID = 1;

% Birth and death rate types
rateSet = {'const-const', '2phase-const'};
rateName = rateSet(rateID);

% Set parameters based on function choices
switch(rateID)
    case 1
        % Constant rate birth-death
        
        % Parameter names and true values
        lam = 0.1;
        mu = 0.5*lam;
        x = [lam mu];
        xdefs = {'lam', 'mu'};
        
        % Points for discretised inference space
        numRV = length(x);
        mi = 30*ones(size(x));
        
        % Space symmetrically centred on true values
        xmax = 10*x;
        xmin = 0.1*x;
        
        % Define rate functions and integrals
        lamt = @(x, tx) x(1)*ones(size(tx));
        mut = @(x, tx) x(2)*ones(size(tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        rhotT = @(x, t1, t2) -(x(1) - x(2))*(t2 - t1);
        
        % Time for which process will run
        T = 80;
        
    case 2
        % Two-phase births and constant death rate

        % Parameter names and true values
        lam1 = 0.01;
        lam2 = 5*lam1;
        mu = 0.5*lam1;
        tau = 70;
        x = [lam1 lam2 tau mu];
        xdefs = {'lam1', 'lam2', 'tau', 'mu'};
        
        % Points for discretised inference space
        numRV = length(x);
        mi = 8*ones(size(x));
        
        % Space symmetrically centred on true values
        xmax = 10*x;
        xmin = 0.1*x;
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(1)*ones(size(tx)).*(tx <= x(3)) + x(2)*ones(size(tx)).*(tx > x(3));
        mut = @(x, tx) x(4)*ones(size(tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        rhotT = @(x, t1, t2) ((x(4) - x(1))*(x(3) - t1) + (x(4) - x(2))*(t2 - x(3)))*(x(3) > t1 && x(3) < t2) +...
            (x(4) - x(1))*(t2 - t1)*(x(3) >= t2) + (x(4) - x(2))*(t2 - t1)*(x(3) <= t1);
        
        % Time for which process will run
        T = 2*tau;
        
    otherwise
        disp('The specified rate id is not supported');
end

% Calculate parameter space
m = prod(mi);
xset = cell(1, 1);
for i = 1:numRV
    xset{i} = linspace(xmin(i), xmax(i), mi(i));
end

% Rate function r(t, s) = e^(int mut-lamt) and integrand of PtT
rst = @(x, t1, t2) exp(rhotT(x, t1, t2));
inttT = @(x, t1, t2) mut(x, t1).*rst(x, t1, t2);


%% Simulate a true birth-death reconstructed tree

% Prob of survival and integral of exp of mut-lamt, conditioned on T
P0T = (1 + integral(@(t2)inttT(x, 0, t2), 0, T))^(-1);
r0T = rst(x, 0, T);

% Calculate prob of no. lineages at T given 1 at 0 (Hohna 2013) and that
% the process survives, Pn0Tsurv
nLinSet = 1:25000;
Pn0Tsurv = zeros(size(nLinSet));
for i = nLinSet
    % From equation (5) and not (3) in Hohna as want process to survive, if
    % used (3) then would have P0T^2 in second term
    Pn0Tsurv(i) = ((1 - P0T*r0T)^(i-1))*(P0T)*r0T;
    %Pn0Tsurv(i) = ((1 - P0T*r0T)^(i-1))*(P0T)*r0T; % equation (3)
end

% Check if simulated enough lineages in nLinSet
if abs(sum(Pn0Tsurv) - 1) < 10^(-10)
    % Normalise and accept
    Pn0Tsurv = Pn0Tsurv/sum(Pn0Tsurv);
else
    % Use larger nLinSet
    disp(['Survival sum is: ' num2str(sum(Pn0Tsurv))]);
    error('Survival probabilities not close enough to 1, increase nLinSet');
end

if ~nBool
    % Generate a uniform r.v. and obtain n from cumsum and nData = n-1
    Pncum = cumsum(Pn0Tsurv);
    n = find(Pncum <= rand, 1, 'last');
    if isempty(n)
        error('n could not be found');
    end
else
    % Set an arbitrary n
    n = 100;
end
disp(['Reconstructed tree has [n T] = ' [num2str(n) ' ' num2str(T)]]);
nLin = 1:n;
nData = n-1;

% Discretise time points for numerical integration
nVals = 20000;
t = linspace(0, T, nVals);

% Calculate PtT from Nee 1994 and the exp rate integral rtT
PtT = zeros(size(t));
rtT = zeros(size(t));
for i = 1:nVals
    PtT(i) = (1 + integral(@(t2)inttT(x, t(i), t2), t(i), T))^(-1);
    rtT(i) = rst(x, t(i), T);
end

% Prob a lineage survives from t to T with exactly 1 descendant (no branch)
p1tT = (PtT.^2).*rtT;

% Draw speciation times for the reconstructed tree using cdf (Hohna 2013)
r = rand(1, n-1);
tden = trapz(t, lamt(x, t).*p1tT, 2);
tnum = cumtrapz(t, lamt(x, t).*p1tT, 2);
tcdf = tnum./tden;

% Use interpolation to find speciation times corresponding to rand values
tspec = interp1(tcdf, t, r);
tspec = sort([0 tspec]);

% Remove duplicate times and set the data length
tspec = unique(tspec);
if length(tspec) ~= n
    error('n distinct speciation times could not be drawn');
end


%% Snyder filtering of the reconstructed process

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
    [tsol, qsol] = ode113(@(ts, y) odeSnyBD(ts, y, xsetMx, numRV, Tsp, nLinCurr, inttT, lamt),...
        [tspec(i) tspec(i+1)], qev(i, :)', options);
    
    % Assign the output values of time and posteriors
    qset{i} = qsol;
    tset{i} = tsol;
    elemLen(i) = length(tsol);
    
    % Perturb the q posterior for the new event
    pert = getNeeRates(tsol(end), xsetMx, numRV, Tsp, nLinCurr, inttT, lamt);
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :).*pert./(qev(i+1, :)*pert');
    
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

%% Analyse and plot results

% Conditional mean
xhat = qn*xsetMx';
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
mu_tn = mut(x, tn);
lam_tn = lamt(x, tn);

% Estimated rates based on function type
switch(rateID)
    case 1
        % Constant rate birth-death (can directly subst cond means)
        muhat = mut(est.xhat, tn);
        lamhat = lamt(est.xhat, tn);

    case 2
        % Two-phase births and constant death rate
        muhat = mut(est.xhat, tn);
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            % Use parameter space and factor in posterior (take mean)
            I1 = t(i) <= xsetMx(3, :);
            lamhat(i) = qnlast*(xsetMx(1, :).*I1 + xsetMx(2, :).*(~I1))';
        end
end

% Collect rates for plotting
rates = [lam_tn lamhat mu_tn muhat];
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


% Plot characteristics of the no. lineages simulated at T
figure;
subplot(2, 1, 1);
plot(nLinSet, Pn0Tsurv);
xlabel('no. lineages at T, n, survival');
ylabel('Pr(N(T) = n | N(0) = 1)');
title(['Probability on no. lineages at T = ' num2str(T)]);
xlim([1 nLinSet(find(Pn0Tsurv < 10^(-6), 1, 'first'))]);
grid;
subplot(2, 1, 2);
stairs(tspec, 1:n);
xlabel('time');
ylabel('no lineages');
title('Lineages through time plot');
grid;

% Plot birth and death rates
figure;
[hAx, hLine1, hLine2] = plotyy(tn, [lam_tn lamhat], tn, [mu_tn muhat]);
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