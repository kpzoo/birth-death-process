% Function to calculate time dependent rate matrix diagonals for the
% time varying birth-death process according to Nee 1994
function ratediag = getNeeTimeVaryInt(t, xsetMx, numRV, Tsp, ncurr)

% Assumptions and modifications
% - uses integral function and Iwasa form
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


% Define net mu-lam rate function from Nee 1994 and its integral
lamb = @(x, tx) x(2)*x(3)./(1 + x(1)*tx);
intB = @(x, t1, t2) x(3)*(t2 - t1) + (x(3)*x(2)/x(1))*(log((x(1)*t1 + 1)./(x(1)*t2 + 1)));

% Constant rate function r(t, s) = e^(int mu-lam) and integrand of PtT
rst = @(t1, t2, x) exp(intB(x, t1, t2));
intT = @(t1, t2, x) lamb(x, t1).*rst(t1, t2, x);

% Calculate P(t, T) from Nee 1994 using exp(r(t, T)) = rst across parameter 
% space for one time set from t to Tsp
PtT = zeros(size(y{1}));
ratediag = zeros(size(y{1}));

% Loop across parameter space and calculate rate diagonal
for k = 1:length(y{1})
    x = [y{1}(k) y{2}(k) y{3}(k)];
    PtT(k) = (1 + integral(@(t2)intT(t, t2, x), t, Tsp))^(-1);
    ratediag(k) = ncurr*lamb(x, t)*PtT(k);
end