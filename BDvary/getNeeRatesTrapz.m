% Function to calculate time dependent rate matrix diagonals for the
% time varying birth-death process according to Nee 1994
function rateDiag = getNeeRatesTrapz(tstart, xsetMx, numRV, Tsp, ncurr, inttT, lamt)

% Assumptions and modifications
% - added inline version of trapz
% - version of getNeeRates which uses trapz
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
rateDiag = zeros(size(y{1}));
tspan = linspace(tstart, Tsp, 200);

% Integrand to be used with trapz to get PtT
%intSet = zeros(size(tspan));

% Loop across parameter space and calculate rate diagonal
for k = 1:length(y{1})
    % Parameter value from space
    xk = xsetMx(1:end, k);
    
    % Get integrand vector for integration
    intSet = inttT(xk, tstart, tspan);
%     for i = 1:length(tspan)
%         intSet(i) = inttT(xk, tstart, tspan(i));
%     end
    % Trapz version of integral for PtT and rate diagonal
    PtT = (1 + trapz(tspan, intSet, 2))^(-1);
    rateDiag(k) = ncurr*lamt(xk, tstart)*PtT;
end

% Remove any negative values of rateDiag (accounts for entries where mut >
% lamt in its parameter space)
rateDiag(rateDiag < 0) = 0;