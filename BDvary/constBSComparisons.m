% Function to perform Paradis and other comparision for constant BD model
function [estParadis, nPts] = constBSComparisons(lam, mu, T, t, tcdf, tspec, xmax, xmin, xdefs, inttT, x)

% Check on limit of PtT as T -> inf in const case
PtTinf = (1 + integral(@(t2)inttT(x, 0, t2), 0, 10^6))^(-1);
PtTinfComp = 1 - mu/lam;
disp(['Check of PtTinf: ' [num2str(PtTinf) ' ' num2str(PtTinfComp)]]);

% Use lots of random draws
rch = rand(1, 10^5);

% Reparameterise the constant rate
sig = lam - mu;
rho = mu/lam;

% Calculate branching times
num = 1 - rho*exp(-sig*T) - rho*(1 - exp(-sig*T))*rch;
den = 1 - rho*exp(-sig*T) - (1 - exp(-sig*T))*rch;
tbranch = (1/sig)*log(num./den);
% Convert times forward from T
tch1 = sort([0 tbranch]);
tch1 = sort(T - tch1);

% Brnching times from Hohna 2013
tch2 = interp1(tcdf, t, rch);
tch2 = sort([0 tch2]);

% Plot comparisons of speciation times
figure;
plot(tch1, tch2);
xlabel('Stadler');
ylabel('Hohner');
title('Comparison of branching times for constant rate birth-death');
grid;

% Calculate Paradis 2011 theoretical cdf
Dfn = @(lam, mu, t, T) exp((lam - mu)*t).*(lam - mu*exp(-(lam - mu)*(T*ones(size(t)) - t))).^(-2);
tDden = trapz(t, Dfn(lam, mu, t, T), 2);
tDnum = cumtrapz(t, Dfn(lam, mu, t, T), 2);
tDcdf = tDnum/tDden;

% Plot comparisons of cdfs of branching times
figure;
plot(tDcdf, tcdf);
xlabel('Paradis');
ylabel('Hohner');
title('Comparison of branching time cdfs for constant rate birth-death');
grid;

% Check T is correct
if T ~= max(tspec)
    error('T needs to be set to max speciation time');
end

% Calculate a Paradis estimate by brute force over same prior range
nPts = 200;
estParadis = testParadis(xmin, xmax, nPts, n, tspec, Dfn, t, T, xdefs, x, tcdf);
disp('Paradis estimate');
disp(estParadis);


% Theoretical constant birth-death speciation distribution - Gernhard2008
fspec = @(l, m, T, y) ((l-m)^2)*(l - m*exp(-(l-m)*T))*exp(-(l-m)*y)...
    /(1 - exp(-(l-m)*T))./((l - m*exp(-(l-m)*y)).^2);

% Plot distribution
yset = 0:0.01:T;
specDistr = fspec(lam, mu, T, yset);

% Check integral is 1
parea = trapz(yset, specDistr);
disp(['Distribution integrates to ' num2str(parea)]);

% Plot distribution
figure;
plot(yset, specDistr);
xlabel('speciation time');
ylabel('probability');
title('Probability density of speciation times for constant birth-death');
grid;

% Fit the distribution to exponential models
ft1 = fittype('a*exp(-b*x)');
ft2 = fittype('a*exp(-a*x)');
[curve1, ~] = fit(yset', specDistr', ft1);
[curve2, ~] = fit(yset', specDistr', ft2);

% Plot fits
gcf;
hold on
plot(curve1, 'c');
plot(curve2);
hold off


%     % Approximate distribution that must be normalised
%     fapprox = @(s, r, y) exp(-s*y)./((1 - r*exp(-s*y)).^2);
%     specApprox = fapprox(mu/lam, lam - mu, yset);
%     specApprox = specApprox/trapz(yset, specApprox);