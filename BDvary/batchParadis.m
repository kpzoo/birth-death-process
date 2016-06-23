% Script to batch compare Paradis and Snyder for constant birth-death
clearvars
clc
close all

tic;

% Define no. replicates and storage variables
M = 5000;
qmarg = cell(1, 1);
xpar = zeros(M, 2);
xsny = zeros(M, 2);

for i = 1:M
    % Obtain Paradis and Snyder estimates
    [est, estParadis, qmarg{i}] = compareParadisSnyConst();
    xpar(i, :) = estParadis.xhat;
    xsny(i, :) = est.xhat;
    disp('/////////////////////////////////////////////////////');
    disp(['Finished batch ' num2str(i) ' of ' num2str(M)]);
    disp('/////////////////////////////////////////////////////');
end
x = est.x;

% Save all data
save('paradisSnyBatch.mat');

% Simulation time
tsim = toc/60;
disp(['Simulation time = ' num2str(tsim) ' mins']);

% Calculate the relative square errors, different parameters on columns
e2par = zeros(size(xpar));
e2sny = zeros(size(xsny));
for i = 1:M
    e2sny(i, :) = 100*(xsny(i, :)./x - 1).^2;
    e2par(i, :) = 100*(xpar(i, :)./x - 1).^2; 
end
Jsny = mean(e2sny);
Jpar = mean(e2par);


% Compare parameter estimates with a boxplot
figure;
subplot(1, 2, 1);
boxplot([xsny(:, 1) xpar(:, 1)], 'labels',{'snyder','paradis'});
ylabel('\lambda');
title(['Estimates of \lambda = ' num2str(x(1))]);
subplot(1, 2, 2);
boxplot([xsny(:, 2) xpar(:, 2)], 'labels',{'snyder','paradis'});
ylabel('\mu');
title(['Estimates of \mu = ' num2str(x(2))]);

% Relative square errors across both estimates
figure;
subplot(2, 2, 1);
hax = histogram(e2sny(:, 1), 20);
hax.FaceAlpha = 0.3; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(Jsny(1)*ones(1, 2), yax, 'r', 'linewidth', 2);
hold off
ylabel('freq');
xlabel('Snyder square errors for \lambda');
legend('square errors', '% MSE', 'location', 'best');
grid;
title(['% MSE for \lambda = ' num2str(Jsny(1))]);

subplot(2, 2, 2);
hax = histogram(e2par(:, 1), 20);
hax.FaceAlpha = 0.3; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(Jpar(1)*ones(1, 2), yax, 'r', 'linewidth', 2);
hold off
ylabel('freq');
xlabel('Paradis square errors for \lambda');
legend('square errors', '% MSE', 'location', 'best');
grid;
title(['% MSE for \lambda = ' num2str(Jpar(1))]);

subplot(2, 2, 3);
hax = histogram(e2sny(:, 2), 20);
hax.FaceAlpha = 0.3; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(Jsny(2)*ones(1, 2), yax, 'r', 'linewidth', 2);
hold off
ylabel('freq');
xlabel('Snyder square errors for \mu');
legend('square errors', '% MSE', 'location', 'best');
grid;
title(['% MSE for \lambda = ' num2str(Jsny(2))]);

subplot(2, 2, 4);
hax = histogram(e2par(:, 2), 20);
hax.FaceAlpha = 0.3; % increase transparency
hold on
% For mean estimate and true value use limits of axes
ax = gca;
yax = ax.YLim;
plot(Jpar(2)*ones(1, 2), yax, 'r', 'linewidth', 2);
hold off
ylabel('freq');
xlabel('Paradis square errors for \mu');
legend('square errors', '% MSE', 'location', 'best');
grid;
title(['% MSE for \mu = ' num2str(Jpar(2))]);


