% Function to find the parameters which minimise the Paradis objective
% function which is also defined here
function paradEst = minParadisNoBoot(xmin, xmax, xset, numRV, n, tspec, tplotcdf, T, lamt,...
    rhotT, inttT, wantPlot, x, xdefs, tcdf)

% Assumptions and modifications
% - removed bootstrap option
% - the empirical cdf matches the matlab function cdfplot(tspec)
% - uses basic constrained minimisation for the cdf, may hit local minima

% Time grid over which integrals are performed and points for internal t(i)
% to T grid
nInt = 2000;
nt = 2000;
t = (0:nt-1).*T/nt;

% Get branching data into form of an empirical cdf (matches cdfplot matlab)
pdata = (1:n)/n;
% Interpolate points based on previous value (zoh) with extrapolation
pinterp = interp1(tspec, pdata, t, 'previous', 'extrap');

% Choose random initial condition from a grid in xset as in Snyder
x0 = zeros(1, numRV);
for i = 1:numRV
    x0(i) = datasample(xset{i}, 1);
end

% Minimisation within bounds of the parameters
[xhat, ssmin] = fmincon(@(xk)paradisObjFn(t, xk, T, nInt, pinterp, lamt, rhotT, inttT),...
    x0, [], [], [], [], xmin, xmax);

% Use estimate to get best cdf from this method
pi = zeros(size(t));
for i = 1:length(t)
    % Define the decreasing time span over which integrals are performed
    tstart = t(i);
    tspan = tstart + (0:nInt-1).*(T - tstart)/nInt;
    
    % Function values with time for integration
    intSet = inttT(xhat, tstart, tspan);
    rho0t = rhotT(xhat, 0, tstart);
    
    % Inline version of trapz for increased speed, assumes row vectors
    trapzIn = diff(tspan)*(intSet(1:nInt-1) + intSet(2:nInt))'/2;
    PtT2 = (1 + trapzIn)^(-2);
    % CDF point value pi(t(i))
    pi(i) = PtT2*exp(-rho0t)*lamt(xhat, tstart);
end
% Get complete CDF for this parameter set
cdf_p = cumtrapz(t, pi, 2)/trapz(t, pi, 2);

% Store estimates and accuracy
paradEst.x = x;
paradEst.xhat = xhat;
paradEst.perc = 100*((1 - xhat./x).^2);
paradEst.lab = xdefs;
paradEst.ssmin = ssmin;
% Store estimated cdf
paradEst.cdf = cdf_p;
paradEst.tcdf = t;

if wantPlot && ~isempty(tcdf)
    % Plot comparisons of CDFs
    figure;
    stairs(tspec, pdata);
    hold on
    plot(tplotcdf, tcdf, 'k');
    plot(t, cdf_p, 'r');
    hold off
    xlabel('time');
    ylabel('cumulative density');
    legend(['emprical, n = ' num2str(n)], 'true CDF',...
        'Paradis CDF', 'location', 'best');
    title('Paradis estimate of time varying birth-death model');
    grid;
end


%% Sub function that performs the constrained minimisation
function ss_err = paradisObjFn(t, xk, T, nInt, pinterp, lamt, rhotT, inttT)

% Calculate Paradis cdf for each point in t
pi = zeros(size(t));
for i = 1:length(t)
    % Define the decreasing time span over which integrals are performed
    tstart = t(i);
    tspan = tstart + (0:nInt-1).*(T - tstart)/nInt;
    
    % Function values with time for integration
    intSet = inttT(xk, tstart, tspan);
    rho0t = rhotT(xk, 0, tstart);
    
    % Inline version of trapz for increased speed, assumes row vectors
    trapzIn = diff(tspan)*(intSet(1:nInt-1) + intSet(2:nInt))'/2;
    PtT2 = (1 + trapzIn)^(-2);
    
    % CDF point value pi(t(i))
    pi(i) = PtT2*exp(-rho0t)*lamt(xk, tstart);
end

% Get complete CDF for this parameter set and calculate sum of square
% errors which is the objective cost
cdf_p = cumtrapz(t, pi, 2)/trapz(t, pi, 2);
ss_err = sum((cdf_p - pinterp).^2);
