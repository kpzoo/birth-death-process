% Function to calculate PtT = = Pr(N(T) > 0 | N(t) = 1) as in Hohna 2013
% for all t in tspan
function PtT = getPtT(lam, mu, tspan)

% Assumptions and modifications
% - modified to include cumtrapz
% - assumes times are synched in tspan
% - inputs are vectors in time of birth and death rates across tspan
% - tspan controls integrals, tspan(1) = t(1), tspan(end) = T = t(end)

% Check input dimensions
if ~(all(size(tspan) == size(lam)) && all(size(tspan) == size(mu)))
    error('Input dimensions inconsistent');
else
    % Ensure row vectors
    if size(lam, 1) ~= 1
        error('Inputs are not row vectors');
    end
end

% Initialise vector
tlen = length(tspan);
PtT = zeros(1, tlen);

% Two loops, inner calculates rtT for changing t, outer does PtT
for i = 1:tlen
    % Partition time interval for integral
    tspace = tspan(i:tlen);
    
    % Calculate all the rtT across this partition
    %rtT = cumtrapz(tspace, mu(i:tlen) - lam(i:tlen), 2);
    
    % Iwasa version
    rtT = cumtrapz(tspace, lam(i:tlen) - mu(i:tlen), 2);
    PtT(i) = (1 + trapz(tspace, lam(i:tlen).*exp(rtT), 2))^(-1);
    
    % Integral with rate difference gives P(t, T), t = tspace(i)
    %PtT(i) = (1 + trapz(tspace, mu(i:tlen).*exp(rtT), 2))^(-1);
end


% Old code which gives same result more slowly
%     rtT = zeros(size(tspace));
%     count = 1;
%     for j = i:tlen
%         % Calculate rate difference r(t, s), t = tspace(1), s = tspace(j)
%         ts = tspan(i:j);
%         rtT(count) = trapz(ts, mu(i:j) - lam(i:j), 2);
%         count = count + 1;
%     end