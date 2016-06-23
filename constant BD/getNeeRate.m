% Function to calculate time dependent rate matrix diagonals for the
% constant birth-death process according to Nee 1994
function lamtdiag = getNeeRate(t, xGen, numRV, Tsp, nLinCurr)

% Assumptions and modifications
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
lamtdiag = nLinCurr*y{2}./(1 - y{1}.*exp(-y{2}*(Tsp - t)));


