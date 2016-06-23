% Script to simulate and estimate time varying birth-death models with
% comparison to analytic results in the constant rate case

% Assumptions and modifications
% - removed data exists boolean and renamed variables
% - does bootstrap version of simBD4 for both Snyder and Paradis
% - includes fmincon version of Paradis estimate
% - can specify the speciation times from a file or generate a new set
% - fixes T to max(tspec) which makes more sense esp for Paradis
% - includes comparison of Paradis cdf estimates for constant rate
% - added again the const rate test of both filter and tree speciations
% - condition on a fixed time, T and n (previously drawed n)
% - uses algorithms from Hohna 2013
% - sampling prob of 1
% - discretised time and used T, mu, a, b (= k in paper) from Nee 1994

clearvars
clc
close all

tic;

% Booleans of importance
useParadis = 1; % to generate a Paradis estimate of parameters
wantBootstrap = 1; % to generate bootstrap estimates for 1 trajectory
wantPlot = 0; % plot from Paradis

%% Initial parameters and definitions

% Define rate function type
rateID = 2;

% Birth and death rate types
rateSet = {'const-const', '2phase-const', 'Nee'};
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
        mi = 20*ones(size(x));
        
        % Space symmetrically centred on true values
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
        mi = 10*ones(size(x));
        
        % Space symmetrically centred on true values
        %xmax = 10*x;
        xmin = 0.1*x;
        xmin(3) = x(3)/2;
        %xmax(3) = x(3)*2;
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(1)*ones(size(tx)).*(tx < x(3)) + x(2)*ones(size(tx)).*(tx >= x(3));
        mut = @(x, tx) x(4)*ones(size(tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        %rhotT = @(x, t1, t2) ((x(4) - x(1))*(x(3) - t1) + (x(4) - x(2))*(t2 - x(3)))*(x(3) > t1 && x(3) < t2) +...
        %   (x(4) - x(1))*(t2 - t1)*(x(3) >= t2) + (x(4) - x(2))*(t2 - t1)*(x(3) <= t1);
        rhotT = @(x, t1, t2) ((x(4) - x(1))*(x(3) - t1) + (x(4) - x(2))*(t2 - x(3)));
        
        % Time for which process will run
        T = 2*tau;
        
        
    case 3
        % Nee 1994 example time varying example
        
        % Set true parameters from Nee 1994
        mu = 0.075;
        k = 4;
        a = 0.2;
        x = [a k mu];
        xdefs = {'a', 'k', 'mu'};
        
        % Space symmetrically centred on true values
        xmin = 0.1*x;
        
        % Points for discretised inference space
        numRV = length(x);
        mi = 20*ones(size(x));
        
        % Define rate functions and integrals, lamt and mut return vectors
        lamt = @(x, tx) x(2)*x(3)./(1 + x(1)*tx);
        mut = @(x, tx) x(3)*ones(size(tx));
        
        % Integral of netRate is rhotT
        netRate = @(x, tx) mut(x, tx) - lamt(x, tx);
        rhotT = @(x, t1, t2) x(3)*(t2 - t1) - (x(2)*x(3)/x(1))*(log(x(1)*t2 + 1) - log(x(1)*t1 + 1));
        
        % Time for which process will run
        T = 100;
        
    otherwise
        disp('The specified rate id is not supported');
end

% Symmetrically place the limits around x
xmax = xmin + 2*x;

% Calculate parameter space
m = prod(mi);
xset = cell(1, 1);
for i = 1:numRV
    xset{i} = linspace(xmin(i), xmax(i), mi(i));
end

% Define no. tips
n = 100;

% Rate function r(t, s) = e^(int mut-lamt) and integrand of PtT
rst = @(x, t1, t2) exp(rhotT(x, t1, t2));
inttT = @(x, t1, t2) mut(x, t2).*rst(x, t1, t2);


%% Simulate a true birth-death reconstructed tree

% Tree tips and coalescent stream length
nLin = 1:n;
nData = n-1;

% Prob of survival and integral of exp of mut-lamt, conditioned on T
P0T = (1 + integral(@(t2)inttT(x, 0, t2), 0, T))^(-1);
r0T = rst(x, 0, T);

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

% Redefine T to maximum speciation time
Torig = T;
T = max(tspec);

% Generate bootstrap data in terms of speciation times sampled from an
% empirical CDF based on tspec
if wantBootstrap
    % No. bootstraps and declared variables
    nBoots = 100;
    tspecB = zeros(nBoots, n); % first entry always set as 0
    
    % Resampling uniform random variables but use 2*n as want to ensure
    % unique times result
    u = rand(nBoots, 2*n);
    u = sort(u, 2);
    
    % Empirical smoothed cdf for interpolation with confidence
    %[pCDF, tspCDF, plb, pub] = ecdf(tspec);
    pCDF = (1:n)/n;
    tspCDF = tspec;
    
    for i = 1:nBoots
        % Find nearest speciation time for cdf values u
        tcheck = unique(sort(interp1(pCDF, tspCDF, u(i, :), 'linear')));
        if length(tcheck) < n
            error(['Not enough unique tspecB values for i = ' num2str(i)]);
        else
            tspecB(i, 2:n) = tcheck(1:n-1);
            tspecB(i, :) = sort(tspecB(i, :));
        end
    end
    % Max speciation times
    TB = max(tspecB, [], 2);
end


%% Paradis estimation

if useParadis
    % Paradis for main data
    paradEst = minParadisNoBoot(xmin, xmax, xset, numRV, n, tspec, t, T, lamt,...
        rhotT, inttT, wantPlot, x, xdefs, tcdf);
    
    if wantBootstrap
        % Variables to store bootstrap estimates
        cdfB = cell(1, nBoots);
        xhatB = zeros(nBoots, numRV);
        percB = zeros(nBoots, numRV);
        ssminB = zeros(1, nBoots);
        
        % Use Paradis method on these bootstrapped times (no distinction for const
        % rate birth-death case)
        for i = 1:nBoots
            paradBoot = minParadisNoBoot(xmin, xmax, xset, numRV, n, tspecB(i, :), t, TB(i), lamt,...
                rhotT, inttT, 0, x, xdefs, []);
            cdfB{i} = paradBoot.cdf;
            xhatB(i, 1:numRV) = paradBoot.xhat;
            percB(i, 1:numRV) = paradBoot.perc;
            ssminB(i) = paradBoot.ssmin;
            disp('**********************************************************');
            disp(['Finished Paradis boot: ' num2str(i) ' of ' num2str(nBoots)]);
            disp('**********************************************************');
        end
        % Simulation time
        tsimPar = toc/60;
        disp(['Simulation time = ' num2str(tsimPar) ' mins']);
        save(['boot_' num2str(rateID) '_' num2str(nBoots)]);
        
        % Histogram of square percentage errrors on each parameter
        mperc = mean(percB, 1);
        figure;
        for i = 1:numRV
            subplot(ceil(numRV/2), 2, i);
            histogram(percB(:, i));
            ylabel('freq')
            xlabel(['% (1 - x_' num2str(i) 'hat/x_' num2str(i) ')^2'])
            title(['Paradis %MSE = ' num2str(mperc(i))]);
            grid;
        end
        
        % Regerate the time points for each cdf
        nt = length(paradBoot.cdf);
        tB = zeros(nBoots, nt);
        for i = 1:nBoots
            tB(i, :) = (0:nt-1).*TB(i)/nt;
        end
        
        % Reduce t to nt values by interpolation
        tred = (0:nt-1).*T/nt;
        tcdfred = interp1(t, tcdf, tred);
        
        % CDF space of bootstrap estimates
        figure;
        plot(tred, tcdfred, 'k', 'linewidth', 2);
        hold on
        h = stairs(tspec, (1:n)/n);
        set(h, 'linewidth', 2);
        plot(paradBoot.tcdf, paradBoot.cdf, 'r', 'linewidth', 2);
        for i = 1:nBoots
            h = plot(tB(i, :), cdfB{i}, ':');
            set(h, 'Color', [0.8 0.8 0.8]);
        end
        hold off
        xlabel('time');
        ylabel('cumulative density');
        legend('true CDF', 'empirical CDF', 'Paradis CDF', 'bootstrap CDFs', 'location', 'best');
        title('Paradis bootstrap for time varying birth-death model');
        grid;
    end
end

%% Snyder filtering of the reconstructed process

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

% Snyder filter original data and obtain conditional mean
[qn, tn] = snyderBootstrap(m, tspec, nData, xsetMx, numRV, rhotT, lamt, mut, nLin);
xhat = qn*xsetMx';

% Display structure for parameters
est.x = x;
est.xhat = xhat(end, :);
est.perc = 100*((1 - est.xhat./x).^2);
est.lab = xdefs;
disp(est);

% Final posterior and marginalisations
qnlast = qn(end, :);
[qmarg, ~] = marginalise(numRV, IDMx, qnlast, mi);

% Bootstrap the Snyder
if wantBootstrap
    % Variables for results
    qnBlast = cell(1, nBoots);
    qmargB = cell(1, nBoots);
    tnB = cell(1, nBoots);
    xhatBS = zeros(nBoots, numRV);
    percBS = zeros(nBoots, numRV);
    % Snyder on each bootstrapped dataset
    for i = 1:nBoots
        [qnB, tnB] = snyderBootstrap(m, tspecB(i, :), nData, xsetMx, numRV, rhotT, lamt, mut, nLin);
        % Last posterior and marginals
        qnBlast{i} = qnB(end, :);
        [qmargB{i}, ~] = marginalise(numRV, IDMx, qnBlast{i}, mi);
        % Conditional estimates
        xhatTemp = qnB*xsetMx';
        xhatBS(i, 1:numRV) = xhatTemp(end, :);
        percBS(i, 1:numRV) = 100*((1 - xhatBS(i, 1:numRV)./x).^2);
        disp('**********************************************************');
        disp(['Finished Snyder boot: ' num2str(i) ' of ' num2str(nBoots)]);
        disp('**********************************************************');
    end
    
    % Simulation time and data store
    tsimSny = toc/60 - tsimPar;
    disp(['Simulation time = ' num2str(tsimSny) ' mins']);
    clear qn qnB
    save(['boot_' num2str(rateID) '_' num2str(nBoots)]);
    
    % Histogram of square percentage errrors on each parameter
    mpercS = mean(percBS, 1);
    figure;
    for i = 1:numRV
        subplot(ceil(numRV/2), 2, i);
        histogram(percBS(:, i));
        ylabel('freq')
        xlabel(['% (1 - x_' num2str(i) 'hat/x_' num2str(i) ')^2'])
        title(['Snyder %MSE = ' num2str(mpercS(i))]);
        grid;
    end
end


%% Analyse and plot results

% Compare boxplots of estimates for each parameter between methods
if wantBootstrap 
    
    % Snyder marginals across bootstraps
    figure;
    for i = 1:numRV
        subplot(ceil(numRV/2), 2, i);
        maxq = 0;
        for j = 1:nBoots
           % Obtain maxima for each posterior
            qval = qmargB{j}{i};
            maxq = max([maxq, max(qval)]);
        end
        % True values
        plot([x(i) x(i)], [0 maxq], 'r', 'linewidth', 2);
        hold on
        for j = 1:nBoots
            % Plot all bootstrap posteriors
            qval = qmargB{j}{i};
            plot(xset{i}, qval, 'b.');
        end
        hold off
        xlabel(['x_' num2str(i)]);
        ylabel(['P(x_' num2str(i) ')']);
        legend('true', 'marginals', 'location', 'best');
        grid;
        title('Bootstrap posteriors');
    end
    
    if useParadis
        
        % Boxplot comparison
        figure;
        for i = 1:numRV
            subplot(1, numRV, i);
            boxplot([xhatBS(:, i) xhatB(:, i)], 'labels',{'snyder','paradis'});
            ylabel(xdefs(i));
            title(['True value = ' num2str(x(i))]);
        end
        
        % Basic comparison
        figure;
        for i = 1:numRV
            subplot(1, numRV, i);
            plot(1:nBoots, [xhatBS(:, i) xhatB(:, i)], 'o');
            hold on
            plot(1:nBoots, x(i)*ones(1, nBoots), 'k');
            hold off
            ylabel(['xhat_' num2str(i)]);
            xlabel('bootstrap no.');
        end
    end

end

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
            I1 = tn(i) <= xsetMx(3, :);
            lamhat(i) = qnlast*(xsetMx(1, :).*I1 + xsetMx(2, :).*(~I1))';
        end
        
    case 3
        % Nee 1994 rates
        muhat = est.xhat(3)*ones(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) = qnlast*(xsetMx(2, :).*xsetMx(3, :)./(1 + xsetMx(1, :)*tn(i)))';
        end
end

% Collect rates for plotting
rates = [lam_tn lamhat mu_tn muhat];
maxr = max(max(rates));
minr = min(min(rates));

% Individual plots for tspec
if wantPlot
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
    
    % Joint final posterior if bivariate
    if numRV == 2
        % Reshape last posterior to 2D form for bivariate plotting
        qtemp = reshape(qnlast, mi);
        if ~isnan(qtemp)
            % Classic joint posterior
            figure;
            meshz(xset{1}, xset{2}, qtemp');
            xlabel('x_1');
            ylabel('x_2');
            zlabel('P(x_1, x_2 | data)');
            title(['Joint posterior with true [x1 x2] = ' num2str(x(1)) ' ' num2str(x(2))]);
            
            % Contour plot
            figure;
            contour(xset{1}, xset{2}, qtemp');
            h = gca;
            xLim = h.XLim;
            yLim = h.YLim;
            hold on
            plot(x(1)*[1 1], yLim, 'k');
            plot(xLim, x(2)*[1 1], 'k');
            xlabel('x_1');
            ylabel('x_2');
            legend('contours', 'true values', 'location', 'best');
            title('Joint posterior contour plot');
            grid;
            
        else
            warning('Mat:qNaN', 'qtemp has NaN values');
        end
    end
    
    
    % Plot characteristics of the no. lineages simulated at T
    figure;
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
end

% Simulation time and data store
tsim = toc/60;
disp(['Simulation time = ' num2str(tsim) ' mins']);
clear qn qnB
save(['bootData_' num2str(rateID) '_' num2str(nBoots)]);