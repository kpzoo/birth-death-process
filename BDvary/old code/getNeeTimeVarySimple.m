% Function to calculate time dependent rate matrix diagonals for the
% time varying birth-death process according to Nee 1994
function ratediag = getNeeTimeVarySimple(t, xsetMx, numRV, Tsp, ncurr)

% Assumptions and modifications
% - only for constant rates at the moment <-------------
% - simplified to use integral function
% - modified to test a constant rate BD process
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
lam1 = y{1};
mu1 = y{2};

% Calculate P(t, T) from Nee 1994 using r(t, T) across parameter space for
% one time set from t to Tsp
PtT = zeros(size(y{1}));
rst = @(tin, uin, mu, lam)exp((mu-lam)*(uin-tin));
intT = @(tin, uin, mu, lam) mu*rst(tin, uin, mu, lam);

% This version from Iwasa works better and matches analytic result
%rst = @(tin, uin, mu, lam)exp((lam-mu)*(uin-tin));
%intT = @(tin, uin, mu, lam) lam*rst(tin, uin, mu, lam);

% Loop across parameter space
for k = 1:length(y{1})
    
    % Iwasa
    %PtT(k) = (1 + integral(@(uin)intT(t, uin, lam1(k), mu1(k)), t, Tsp))^(-1);
    %PtT(k) = (1 + integral(@(tin)intT(tin, Tsp, lam1(k), mu1(k)), t, Tsp))^(-1);
    
    % Nee
    PtT(k) = (1 + integral(@(uin)intT(t, uin, lam1(k), mu1(k)), t, Tsp))^(-1);
end

% Calculate reconstructed conditioned birth rate from Nee 1994
ratediag = ncurr*lam1.*PtT;