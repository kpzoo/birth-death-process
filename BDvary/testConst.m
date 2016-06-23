% Function to analytically Snyder filter constant rate birth-death
function [est2, qmarg2] = testConst(x, tspec, m, nLin, nData, xset, q0,...
    xsetMx, numRV, Tsp, xdefs, IDMx, mi)


% Posterior vectors on events
qev2 = zeros(nData+1, m);
qev2(1, :) = q0;
tev2 = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset2 = cell(1, 1);
tset2 = cell(1, 1);
options = odeset('NonNegative', 1:m);
elemLen2 = zeros(1, nData);

% Loop across coalescent events and perform filtering 
for i = 1:nData
    % Current no. lineages
    nLinCurr = nLin(i);
    
    % Solve linear ODEs continuously with setting of options, no Q
    [tsol2, qsol2] = ode113(@(ts, y) odeSnyBDconstLamMu(ts, y, xsetMx, numRV, Tsp, nLinCurr),...
        [tspec(i) tspec(i+1)], qev2(i, :)', options);
    
    % Assign the output values of time and posteriors
    qset2{i} = qsol2;
    tset2{i} = tsol2;
    elemLen2(i) = length(tsol2);
    
    % Perturb the q posterior for the new event
    pert2 = getNeeRateLamMu(tsol2(end), xsetMx, numRV, Tsp, nLinCurr);
    qev2(i+1, :) = qsol2(end, :);
    tev2(i+1) = tsol2(end);
    qev2(i+1, :) = qev2(i+1, :).*pert2./(qev2(i+1, :)*pert2');
    
    disp(['Finished: ' num2str(i) ' of ' num2str(nData)]);
    
end

% Get full length of ODE solution data and assign appending vectors
lenFull2 = sum(elemLen2);
stop = 0;
qn2 = -ones(lenFull2, m);
tn2 = -ones(lenFull2, 1);

% Append the cell based posterior and time data into a single structure
for i = 1:nData
    % Loop calculates the start and end points along the array successively
    % and then assigns the appropriate cell element
    start = stop + 1;
    stop = elemLen2(i) + start - 1;
    tn2(start:stop) = tset2{i};
    qn2(start:stop, :) = qset2{i};
end

% Conditional mean
xhat2 = qn2*xsetMx';

% Display structure for parameters
est2.x = x;
est2.xhat = xhat2(end, :);
est2.perc = 100*(1 - est2.xhat./x);
est2.lab = xdefs;
disp(est2);

% Final posterior and marginalisations
qnlast2 = qn2(end, :);
[qmarg2, ~] = marginalise(numRV, IDMx, qnlast2, mi);

% Plot marginals for each parameter and real value
figure;
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    plot(xset{i}, qmarg2{i}, 'b', [x(i) x(i)], [0 max(qmarg2{i})], 'k',...
        [xhat2(end, i) xhat2(end, i)], [0 max(qmarg2{i})], 'r');
    xlabel(['x_' num2str(i)]);
    ylabel(['P(x_' num2str(i) ')']);
    legend('marginal', 'true', 'cond mean', 'location', 'best');
    grid;
end
