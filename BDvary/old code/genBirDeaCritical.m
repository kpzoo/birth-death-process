% Function to sample from a critical birth-death process according to
% Hartmann et al 2008 where lam = mu
function [T, tbranch] = genBirDeaCritical(n, mu)

% Assumptions and modifications
% - inputs - no. extant lineages, n
% - outputs - speciation times
% - only for critical birth-death process

% Sample n uniform variables in [0, 1] for (n-1) times and one origin
r = rand(1, n-1);
r0 = rand;

% Draw a tree age
T = 1/(mu*(r0^(-1/n) - 1));

% Calculate branching times
tbranch = r*T./(1 + mu*T*(1 - r));
