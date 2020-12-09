function [gamma] = util_FitExponential(Q, t, fitting_type)
%util_FitExponential fits an exponential function to recession segments.
%   a: Q = a+b*exp(-gamma*t) (stable when using lsqnonlin)
%   b: Q = Q0*exp(-gamma*t) (seems to be the most stable when using fitnlm)
%
%   INPUT
%   Q: dependent variable, typically streamflow in [mm/timestep]
%   t: independent variable, typically time
%   OPTIONAL
%   fitting_type: what kind of exponential function to fit ('a' or 'b').
%
%   OUTPUT
%   gamma: fitted time parameter of exponential function
%
%   EXAMPLE
%   t = [0:10]';
%   Q = exp(-0.1*t);
%   [gamma] = util_FitExponential(Q, t)
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 3
    fitting_type = 'a';
end

switch fitting_type
    case 'a'
        t = [0:length(t)-1]'; % start of recession equals start of exponential
        % lsqnonlin is more stable than fitnlm but requires optimization toolbox
        ExponentialObjective = @(para) para(2) + para(3).*exp(-para(1).*t) - Q;
        para0 = [0.1 0.1 0.1];
        options = optimoptions(@lsqnonlin,'Display','off');
        para = lsqnonlin(ExponentialObjective, para0, [1e-9 1e-9 1e-9], [100 100 100], options);
        gamma = para(1);
        
    case 'b'
        t = [0:length(t)-1]'; % start of recession equals start of exponential
        % lsqnonlin is more stable than fitnlm but requires optimization toolbox
        ExponentialObjective = @(para) Q(1).*exp(-para(1).*t) - Q;
        para0 = [0.1];
        options = optimoptions(@lsqnonlin,'Display','off');
        para = lsqnonlin(ExponentialObjective, para0, [1e-9], [100], options);
        gamma = para(1);
        
    otherwise
        error('Please choose one of the available fitting types: a or b.')
end

%{
% plot
figure; plot(t,Q,'o'); hold on
% Q_est = Q(1)*exp(t.*-gamma);
a = para(2);
b = para(3);
gamma = para(1);
Q_est = a+b*exp(t.*-gamma);
plot(t,Q_est)
%}

end