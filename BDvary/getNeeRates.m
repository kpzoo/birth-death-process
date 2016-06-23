% Function to calculate time dependent rate matrix diagonals for the
% time varying birth-death process according to Nee 1994
function rateDiag = getNeeRates(tstart, xsetMx, numRV, Tsp, ncurr, inttT, lamt)

% Assumptions and modifications
% - uses function handle input
% - simplified to use integral function
% - Tsp is max speciation time
% - time input must be a scalar

% Get variables in form for direction function input
y = cell(1, numRV);
for i = 1:numRV
        % Want space of rate so use xsetMx
        y{i} = xsetMx(i, :);
end

% Calculate P(t, T) from Nee 1994 using exp(r(t, T)) = rst across parameter 
% space for one time set from t to Tsp
PtT = zeros(size(y{1}));
rateDiag = zeros(size(y{1}));

% Loop across parameter space and calculate rate diagonal
for k = 1:length(y{1})
    xk = xsetMx(1:end, k);
    PtT(k) = (1 + integral(@(t2)inttT(xk, tstart, t2), tstart, Tsp))^(-1);
    rateDiag(k) = ncurr*lamt(xk, tstart)*PtT(k);
end

% Remove any negative values of rateDiag (accounts for entries where mut >
% lamt in its parameter space)
rateDiag(rateDiag < 0) = 0;