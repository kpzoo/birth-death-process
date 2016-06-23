% Script to batch testBD for constant rate birth-death process

% Assumptions and modifications
% - uses functional form of testBD


clc
close all
clearvars

% Set simulation parameters
M = 10000;
n = 200;
mi = [100 100];
numRV = length(mi);

% Set space for estimated parameters
xdefs = {'rho', 'sigma'};
%minSpace = [0.01 0.01]; % original lam and mu
%maxSpace = [0.99 100];

minSpace = [0.01 0.01]; % low mu case
maxSpace = [0.1 2];

% Get the instance of the random variables from the parameter space
xset = cell(1, 1);
x = zeros(1, numRV);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
    % Randomly sample uniformly across range
    %x(i) = datasample(xset{i}, 1);
    % Sample at mid point of range
    x(i) = xset{i}(round(mi(i)/2));
end

% True parameters with names
rho = x(1);
sig = x(2);
lam = sig/(1 - rho);
mu = lam - sig;

% Storage parameters
rho_h = zeros(1, M);
sig_h = rho_h;
lam_h = rho_h;
mu_h = rho_h;

% Run M iterations with key parameter taken at mid point
for i = 1:M
    [birdea, est] = testBDfn(n, mi, xdefs, xset, x);
    % Parameter estimates
    rho_h(i) = est.xhat(1);
    sig_h(i) = est.xhat(2);
    lam_h(i) = birdea.est(1);
    mu_h(i) = birdea.est(2);
    disp(['Finished ' num2str(i) ' of ' num2str(M)]);
end

% Save data
save(['batBDlowMu_' num2str(n) '_' num2str(M)], 'rho', 'sig', 'lam', 'mu',...
    'rho_h', 'sig_h', 'lam_h', 'mu_h', 'M', 'n', 'mi');

% Raw errors
rho_e = rho - rho_h;
sig_e = sig - sig_h;
lam_e = lam - lam_h;
mu_e = mu - mu_h;

% Mean estimates
rho_av = mean(rho_h);
sig_av = mean(sig_h);
lam_av = mean(lam_h);
mu_av = mean(mu_h);

% Mean square errors
rho_sq = rho_e.^2;
rho_mse = mean(rho_sq);
sig_sq = sig_e.^2;
sig_mse = mean(sig_sq);
lam_sq = lam_e.^2;
lam_mse = mean(lam_sq);
mu_sq = mu_e.^2;
mu_mse = mean(mu_sq);

% Mean square percentage errors
rho_perc = 100*(1 - rho_h/rho).^2;
rho_m = mean(rho_perc);
sig_perc = 100*(1 - sig_h/sig).^2;
sig_m = mean(sig_perc);
lam_perc = 100*(1 - lam_h/lam).^2;
lam_m = mean(lam_perc);
mu_perc = 100*(1 - mu_h/mu).^2;
mu_m = mean(mu_perc);

% Display mean percentage errors
perc.rho = rho_m;
perc.sig = sig_m;
perc.lam = lam_m;
perc.mu = mu_m;
disp(perc);

% Plot histograms with mean and true value
figure;
% Parameter: rho
subplot(1, 2, 1);
hax = histogram(rho_h);
hax.FaceAlpha = 0.03; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(rho_av*ones(1, 2), yax, 'r', rho*ones(1, 2), yax, 'k', 'linewidth', 3);
hold off
ylabel('freq');
xlabel('\rho');
legend('estimates', 'mean est', 'true', 'location', 'best');
grid;
title(['% MSE = ' num2str(rho_m)]);
% Parameter: sig
subplot(1, 2, 2);
hax = histogram(sig_h);
hax.FaceAlpha = 0.03; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(sig_av*ones(1, 2), yax, 'r', sig*ones(1, 2), yax, 'k', 'linewidth', 3);
hold off
ylabel('freq');
xlabel('\sigma');
legend('estimates', 'mean est', 'true', 'location', 'best');
grid;
title(['% MSE = ' num2str(sig_m)]);

% Plot histograms with mean and true value
figure;
% Parameter: lam
subplot(1, 2, 1);
hax = histogram(lam_h);
hax.FaceAlpha = 0.03; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(lam_av*ones(1, 2), yax, 'r', lam*ones(1, 2), yax, 'k', 'linewidth', 3);
hold off
ylabel('freq');
xlabel('\lambda');
legend('estimates', 'mean est', 'true', 'location', 'best');
grid;
title(['% MSE = ' num2str(lam_m)]);
% Parameter: mu
subplot(1, 2, 2);
hax = histogram(mu_h);
hax.FaceAlpha = 0.03; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(mu_av*ones(1, 2), yax, 'r', mu*ones(1, 2), yax, 'k', 'linewidth', 3);
hold off
ylabel('freq');
xlabel('\mu');
legend('estimates', 'mean est', 'true', 'location', 'best');
grid;
title(['% MSE = ' num2str(mu_m)]);