% Function to calculate time varying Paradis 2011 estimate
function estParadis = tvaryParadis(xmin, xmax, nPtperVar, numRV, n, tspec, tplotcdf, T, lamt,...
    rhotT, inttT, wantPlot, x, xdefs, tcdf)

% Assumptions and modifications
% - added inline functions for linspace and trapz
% - need correction for lamt - mut < 0?
% - internal t spanning 0 to T
% - uses trapz integrals
% - a fine m has to be used since it defines the resolution of the method
% - this m and mi are not the same as in Snyder but xset is identical

%profile on

% Check T is correct
if T ~= max(tspec)
    error('T must be the largest speciation time');
end

% Time span that defines the integrals of trapz
% t = linspace(0, T, 1000);
% Inline linspace
t = (0:999).*T/1000;

% Get branching data into form of an empirical cdf
pdata = (1:n)/n;
% Interpolate points based on previous value (zoh) with extrapolation
pinterp = interp1(tspec, pdata, t, 'previous', 'extrap');

% Calculate a finer parametric search space of xsetMx form
mi = nPtperVar*ones(1, numRV);
m = prod(mi);
xset = cell(1, 1);
for i = 1:numRV
    xset{i} = linspace(xmin(i), xmax(i), mi(i));
end

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

% Storage of square errors and cdfs
ss_err = zeros(1, m);
cdf_p = cell(1, 1);
nInt = 100; % for integration

% Loop across all possible parameter sets and
for k = 1:m
    % Current parameter set being evaluated
    xk = xsetMx(1:end, k);
    
    % Calculate Paradis cdf for each point in t
    pi = zeros(size(t));
    for i = 1:length(t)
        tstart = t(i);
        %tspan = linspace(tstart, T, 100);
        tspan = tstart + (0:nInt-1).*(T - tstart)/nInt;
        intSet = inttT(xk, tstart, tspan);
        rho0t = rhotT(xk, 0, tstart);
        %PtT2 = (1 + trapz(tspan, intSet, 2))^(-2);
        
        % Inline version of trapz for increased speed, assumes row vectors
        trapzIn = diff(tspan)*(intSet(1:nInt-1) + intSet(2:nInt))'/2;
        PtT2 = (1 + trapzIn)^(-2);
        % CDF point value pi(t(i))
        pi(i) = PtT2*exp(-rho0t)*lamt(xk, tstart);
    end
    
    % Get complete CDF for this parameter set
    cdf_p{k} = cumtrapz(t, pi, 2)/trapz(t, pi, 2);
    ss_err(k) = sum((cdf_p{k} - pinterp).^2);
    disp(['Finished ' num2str(k) ' of ' num2str(m)]);
end

% Get arguments which minimise ss_err => estimates
idest = find(ss_err == min(min(ss_err)));
% Store estimates from Paradis
estParadis.x = x;
estParadis.xhat = xsetMx(1:end, idest)';
estParadis.perc = 100*(1 - estParadis.xhat./x);
estParadis.lab = xdefs;
estParadis.ss_err = ss_err;
estParadis.cdf = cdf_p{idest};
estParadis.xsetMx = xsetMx;

if wantPlot
    % Plot comparisons
    figure;
    stairs(tspec, pdata);
    hold on
    plot(tplotcdf, tcdf, 'k');
    plot(t, cdf_p{idest}, 'r');
    xlabel('time');
    ylabel('cumulative density');
    legend(['emprical, n = ' num2str(n)], 'true CDF',...
        'Paradis CDF', 'location', 'best');
    title('Paradis estimate of time varying birth-death model');
    grid;
    
    % Plot sum of squares across m combinations
    figure;
    plot(1:m, ss_err);
    xlabel('points');
    ylabel('sum of square errors');
    
end

%profile viewer
