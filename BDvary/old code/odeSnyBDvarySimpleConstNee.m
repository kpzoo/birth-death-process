% Function to setup Snyder ODE set for constant rate birth-death

% Assumptions
% - uses proper Nee
% - uses simplified integral case which assumes constant
% - modified to test constant rate case
% - assumes time varying birth-death from Nee 1994
% - uses equivalent birth rate for rate matrix from Nee 1994
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP
% - binfac is the appropriate binomial factor based on no. events

function dy = odeSnyBDvarySimpleConstNee(ts, y, xsetMx, numRV, Tsp, nLinCurr)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Diagonal of time dependent rate matrix for reconstructed birth (Nee 1994)
lamtdiag = getNeeTimeVarySimpleNee(ts, xsetMx, numRV, Tsp, nLinCurr);

% Solve non-linear differential equation set - RV filtering
nonLinDiag = lamtdiag*y;
dy = y'.*(-lamtdiag + nonLinDiag);

% Ensure output is column vector assuming input was
dy = dy';