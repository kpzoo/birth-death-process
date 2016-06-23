% Function to calculate time dependent rate matrix diagonals for the
% time varying birth-death process according to Nee 1994
function ratediag = getNeeTimeVaryIntFn(t, xsetMx, numRV, Tsp, ncurr, intT, lamb)

% Assumptions and modifications
% - uses function handle input of mub instead of lamb in PtT
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

% Calculate P(t, T) from Nee 1994 using exp(r(t, T)) = rst across parameter 
% space for one time set from t to Tsp
PtT = zeros(size(y{1}));
ratediag = zeros(size(y{1}));

% Loop across parameter space and calculate rate diagonal
for k = 1:length(y{1})
    x = xsetMx(1:end, k);
    PtT(k) = (1 + integral(@(t2)intT(t, t2, x), t, Tsp))^(-1);
    ratediag(k) = ncurr*lamb(x, t)*PtT(k);
end