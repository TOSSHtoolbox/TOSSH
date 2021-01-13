function [gamma] = util_FitExponential(Q, t, fitting_type)
%util_FitExponential fits an exponential function to recession segments.
%   Different types of functions and fitting options are available.
%   Q = Q0*exp(-gamma*t) (either linear regression in semilog space or
%       nonlinear fit)
%   Q = b*exp(-gamma*t) (nonlinear two-parameter fit)
%   Q = a+b*exp(-gamma*t) (nonlinear three-parameter fit)
%   Note that the non-linear fitting methods are much slower.
%
%   INPUT
%   Q: dependent variable, typically streamflow in [mm/timestep]
%   t: independent variable, typically time
%   OPTIONAL
%   fitting_type: what kind of exponential function to fit ('semilog',
%       'nonlinear', 'nonlinear2', 'nonlinear3')
%
%   OUTPUT
%   gamma: fitted parameter of exponential function, typically [1/timestep]
%
%   EXAMPLE
%   t = [0:10]';
%   Q = exp(-0.1*t);
%   [gamma] = util_FitExponential(Q, t);
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 3
    fitting_type = 'semilog';
end

t = [0:length(t)-1]'; % start of recession equals start of exponential

switch fitting_type
    
    case 'semilog'
        Q = Q(:); % make sure that Q is a column vector
        gamma = -t\(log(Q)-log(Q(1)));
        
    case 'nonlinear'
        ExponentialObjective = @(para) Q(1).*exp(-para(1).*t) - Q;
        para0 = [0.1];
        options = optimoptions(@lsqnonlin,'Display','off');
        para = lsqnonlin(ExponentialObjective, para0, [1e-6], [100], options);
        gamma = para(1);
        
    case 'nonlinear2'
        ExponentialObjective = @(para) para(2).*exp(-para(1).*t) - Q;
        para0 = [0.1 0.1];
        options = optimoptions(@lsqnonlin,'Display','off');
        para = lsqnonlin(ExponentialObjective, para0, [1e-6 1e-6], [100 100], options);
        gamma = para(1);
        
    case 'nonlinear3'
        ExponentialObjective = @(para) para(2) + para(3).*exp(-para(1).*t) - Q;
        para0 = [0.1 0.1 0.1];
        options = optimoptions(@lsqnonlin,'Display','off');
        para = lsqnonlin(ExponentialObjective, para0, [1e-6 1e-6 1e-6], [100 100 100], options);
        gamma = para(1);
        
    otherwise
        error('Please choose one of the available fitting types: a or b.')
end


% plot
%{
figure; plot(t,Q,'o'); hold on
switch fitting_type
    case 'nonlinear2'
        b = para(2);
        gamma = para(1);
        Q_est = b*exp(t.*-gamma);
    case 'nonlinear3'
        a = para(2);
        b = para(3);
        gamma = para(1);
        Q_est = a+b*exp(t.*-gamma);
    otherwise
        Q_est = Q(1)*exp(t.*-gamma);
end
plot(t,Q_est)
%}

end