% Function to save the variables that are important for other simulations
function data = saveSimData(name, est, mi, n , Torig, Tsp, T, tspec, xmax,...
    xmin, xsetMx, IDMx, t, tcdf, xset, x)


% Important variables that define the simulation and results
data.est= est;
data.ics = [n Torig Tsp T];
data.icsLab = {'n', 'Torig', 'Tsp', 'T'};
data.mi = mi;
data.tspec = tspec;
data.xmax = xmax;
data.xmin = xmin;
data.xsetMx = xsetMx;
data.IDMx = IDMx;
data.t = t;
data.tcdf = tcdf;
data.xset = xset;
data.x = x;

save(name, 'data');
