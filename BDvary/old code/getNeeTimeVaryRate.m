% Function to calculate time dependent rate matrix diagonals for the
% time varying birth-death process according to Nee 1994
function ratediag = getNeeTimeVaryRate(t, xsetMx, numRV, Tsp, ncurr)

% Assumptions and modifications
% - Nee example time varying birth death process
% - Tsp is max speciation time
% - time input must be a scalar

% Get variables in form for direction function input
y = cell(1, numRV);
for i = 1:numRV
        % Want space of rate so use xsetMx
        y{i} = xsetMx(i, :);
end


% Birth and death rate set with Nee 1994 rates
lam1 = y{2}.*y{3}./(1 + t*y{1});
mu1 = y{3};

% Calculate P(t, T) from Nee 1994 using r(t, T) across parameter space
tspace = linspace(t, Tsp, 1000);
%rtT = zeros(size(tspace));
PtT = zeros(size(y{1}));
% Outer loop across parameter space, inner across time for rtT
for k = 1:length(y{1})
    % Calculate all the rtT across this partition
    %rtT = cumtrapz(tspace, mu1(k) - y{2}(k)*y{3}(k)./(1 + tspace*y{1}(k)), 2);
    
    % Iwasa
    rtT = cumtrapz(tspace, -(mu1(k) - y{2}(k)*y{3}(k)./(1 + tspace*y{1}(k))), 2);
    
    % Integral with rate difference gives P(t, T), t = tspace(i)
    %PtT(k) = (1 + trapz(tspace, mu1(k)*exp(rtT), 2))^(-1);
    
    % Iwasa
    PtT(k) = (1 + trapz(tspace, lam1(k)*exp(rtT), 2))^(-1);
    
end


% Calculate reconstructed conditioned birth rate from Nee 1994
ratediag = ncurr*lam1.*PtT;


% Old code
% for j = 1:100
%     % Calculate rate difference r(t, s), t = tspace(1), s = tspace(j)
%     ts = tspace(1:j);
%     % Calculate birth and death rates across ts for each parameter val
%     dr = mu1(k) - y{2}(k)*y{3}(k)./(1 + ts*y{1}(k));
%     rtT(j) = trapz(ts, dr, 2);
% end