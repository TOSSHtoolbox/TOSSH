function [a, b, error_flag, error_str] = ...
    util_FitPowerLaw(Q_rec, dQdt_rec, fitting_type, weights)
%util_FitPowerLaw fits a power law (to recession segments).
%   dQ/dt = -a*Q^b
%   Different options are available (e.g. linear fitting in loglog-space,
%   non-linear, or fitting a fixed slope).
%
%   INPUT
%   Q_rec: streamflow (recession periods only)
%   dQdt_rec: corresponding flow rate gradients
%   fitting_type: specifies fitting procedure: linear, non-linear, linear
%       with slope 1, linear with slope 2
%   weights: weights to fit curve
%
%   OUTPUT
%   a: scaling parameter
%   b: parameter of non-linearity
%   error_flag: 0 (no error), 1 (warning), 2 (errorin data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   flow_section = util_RecessionSegments(Q,t,'plot_results',true); % get recession segments
%   [dQdt, Qm, flow_section, R2] = util_dQdt(Q, t, flow_section); % get flow rate gradient
%   rec = ~isnan(Qm);
%   [a, b] = util_FitPowerLaw(Qm(rec), dQdt(rec));
%   [a, b] = util_FitPowerLaw(Qm(rec), dQdt(rec), 'fitting_type', 'linear', 'weights', R2(rec))
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 2
    error('Not enough input arguments.')
elseif nargin < 3
    fitting_type = 'linear';
elseif nargin < 4
    weights = ones(size(Q_rec));
end

% default setting reads as good data
error_flag = 0;
error_str = '';

% data checks
if all(isnan(Q_rec)) || all(isnan(dQdt_rec))
    a = NaN;
    b = NaN;
    error_flag = 1;
    error_str = ['Warning: Some recessions consist only of NaN values (possibly because eps > 0). ', error_str];
    return
end

if strcmp(fitting_type, 'linear') % weighted linear regression in log log space
    
    A = [ones(size(Q_rec)), log(Q_rec)];
    P = (weights.*A)\(weights.*log(-dQdt_rec));
    a = exp(P(1));
    b = P(2);
    
    %{
    linFcn = @(p,x) p(2).*x + p(1);
    p0 = [0.1, 1.0];
    fit_nonlin = fitnlm(log(Q_rec),log(-dQdt_rec),linFcn,p0,'Weight',weights);
    a = exp(fit_nonlin.Coefficients.Estimate(1));
    b = fit_nonlin.Coefficients.Estimate(2);
    %}
    
elseif strcmp(fitting_type, 'nonlinear') % weighted nonlinear regression in lin space
    
    powerFcn = @(p,x) p(1).*x.^p(2);
    p0 = [0.1, 1.0];
    fit_nonlin = fitnlm(Q_rec,-dQdt_rec,powerFcn,p0,'Weight',weights);
    a = fit_nonlin.Coefficients.Estimate(1);
    b = fit_nonlin.Coefficients.Estimate(2);
    
elseif strcmp(fitting_type, 'slope1') % linear regression in log log space with fixed slope 1
    
    A = [ones(size(Q_rec))];
    b = 1;
    B = log(-dQdt_rec) - b*log(Q_rec);
    a = exp((weights.*A)\(weights.*B));
    
elseif strcmp(fitting_type, 'slope2') % linear regression in log log space with fixed slope 2
    
    A = [ones(size(Q_rec))];
    b = 2;
    B = log(-dQdt_rec) - b*log(Q_rec);
    a = exp((weights.*A)\(weights.*B));
    
else
    error('Invalid fitting type.')
end

end


