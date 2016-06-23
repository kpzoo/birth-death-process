% Function to obtain the Paradis 2011 estimate for constant birth-death
function estParadis = testParadis(xmin, xmax, nPts, n, tspec, Dfn, t, T, xdefs, x, tcdf)


    % Use the same space limits but on a finer scale from nPts
    lam_p = linspace(xmin(1), xmax(1), nPts);
    mu_p = linspace(xmin(2), xmax(2), nPts);
    
    % Get branching data into form of an empirical cdf
    pdata = (1:n)/n;
    % Interpolate points based on previous value (zoh) with extrapolation
    pinterp = interp1(tspec, pdata, t, 'previous', 'extrap');
    
    % Calculate sum of square errors
    ss_err = zeros(nPts, nPts);
    for i = 1:nPts
        for j = 1:nPts
            if lam_p(i) > mu_p(j)
                den_p = trapz(t, Dfn(lam_p(i), mu_p(j), t, T), 2);
                num_p = cumtrapz(t, Dfn(lam_p(i), mu_p(j), t, T), 2);
                cdf_p = num_p/den_p;
                ss_err(i, j) = sum((cdf_p - pinterp).^2);
            else
                ss_err(i, j) = inf;
            end
        end
        disp(['Finished Paradis run ' num2str(i) ' of ' num2str(nPts)]);
    end
    
    % Get arguments which minimise ss_err => estimates
    [ival, jval] = find(ss_err == min(min(ss_err)));
    lam_ph = lam_p(ival);
    mu_ph = mu_p(jval);
    
    % Store estimates from Paradis
    estParadis.x = x;
    estParadis.xhat = [lam_ph mu_ph];
    estParadis.perc = 100*(1 - estParadis.xhat./x);
    estParadis.lab = xdefs;

%     % Calculate the best estimate's cdf
%     den_p = trapz(t, Dfn(lam_ph, mu_ph, t, T), 2);
%     num_p = cumtrapz(t, Dfn(lam_ph, mu_ph, t, T), 2);
%     tcdfParadis = num_p/den_p;
%             
%     % Plot comparisons
%     figure;
%     stairs(tspec, pdata);
%     hold on
%     plot(t, tcdf, 'k');
%     plot(t, tcdfParadis, 'r');
%     xlabel('time');
%     ylabel('cumulative density');
%     legend(['emprical, n = ' num2str(n)], ['true, [\lambda \mu] = ' [num2str(x(1)) ' ' num2str(x(2))]],...
%         ['Paradis est, [\lambda \mu] = ' [num2str(lam_ph) ' ' num2str(mu_ph)]], 'location', 'best');
%     title('Paradis estimate of constant rate birth-death model');
%     grid;