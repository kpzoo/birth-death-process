% Function to calculate time dependent rate matrix diagonals for the
% constant birth-death process according to Nee 1994
function lamtdiag = getNeeRateLamMu(t, xGen, numRV, Tsp, nLinCurr)

% Assumptions and modifications
% - version of getNeeRate directly for lam and mu as y{1} and y{2}
% - Tsp is max speciation time
% - time input must be a scalar


% Get variables in form for direction function input
y = cell(1, numRV);
for i = 1:numRV
        % Want space of lamt so use xsetMx or column form of x if want true
        % value instead of estimates
        y{i} = xGen(i, :);
end

% Calculate reconstructed conditioned birth rate from Nee 1994
%sigma = y{1} - y{2};
%rho = y{2}./y{1};
%lamtdiag = nLinCurr*sigma./(1 - rho.*exp(-sigma*(Tsp - t)));


lamtdiag = nLinCurr*(y{1} - y{2})./(1 - (y{2}./y{1}).*exp(-(y{1} - y{2})*(Tsp - t)));


