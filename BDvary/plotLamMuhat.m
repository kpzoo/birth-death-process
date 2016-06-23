% Function for plotting estimated rates
function [lamhat, muhat] = plotLamMuhat(xsetMx, xhat, tn, qnlast, rateID, mut, lamt)

% Obtain estimates of parametric birth/death reates from posterior
switch(rateID)
    case 1
        % Constant rate birth-death (can directly subst cond means)
        muhat = mut(xhat, tn);
        lamhat = lamt(xhat, tn);
        
    case 2
        % Two-phase births and constant death rate
        muhat = mut(xhat, tn);
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            % Use parameter space and factor in posterior (take mean)
            I1 = tn(i) <= xsetMx(3, :);
            lamhat(i) = qnlast*(xsetMx(1, :).*I1 + xsetMx(2, :).*(~I1))';
        end
        
    case 3
        % Nee 1994 rates
        muhat = xhat(3)*ones(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) = qnlast*(xsetMx(2, :).*xsetMx(3, :)./(1 + xsetMx(1, :)*tn(i)))';
        end
        
    case 4
        % Paradis logistic rates
        muhat = zeros(size(tn));
        lamhat = zeros(size(tn));
        for i = 1:length(tn)
            lamhat(i) = qnlast*(1./(1 + exp(-xsetMx(1, :)*tn(i) + xsetMx(2, :))))';
            muhat(i) = qnlast*(1./(1 + exp(-xsetMx(3, :)*tn(i) + xsetMx(4, :))))';
        end
end