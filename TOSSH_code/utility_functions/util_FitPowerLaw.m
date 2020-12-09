function [a, b, error_flag, error_str] = util_FitPowerLaw(Q_rec, dQdt_rec, varargin)
%util_FitPowerLaw fits a power law (to recession segments).
%   dQ/dt = - a Q ^ b
%   Different options are available (e.g. linear fitting in loglog-space,
%   non-linear, or fitting a fixed slope).
%
%   INPUT
%   Q_rec: streamflow (recession periods only)
%   dQdt_rec: corresponding flow rate gradients
%   OPTIONAL
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
end

ip = inputParser;
ip.CaseSensitive = true;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q_rec', @(Q_rec) isnumeric(Q_rec) && (size(Q_rec,1)==1 || size(Q_rec,2)==1))
addRequired(ip, 'dQdt_rec', @(dQdt_rec) isnumeric(dQdt_rec) && (size(dQdt_rec,1)==1 || size(dQdt_rec,2)==1)) 

% optional input arguments
addParameter(ip, 'fitting_type', 'linear', @ischar) %
addParameter(ip, 'weights', ones(size(Q_rec)), @isnumeric) % 

parse(ip, Q_rec, dQdt_rec, varargin{:})
fitting_type = ip.Results.fitting_type;
weights = ip.Results.weights;

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

if strcmp(fitting_type, 'linear') % linear regression in log log space
    
    linFcn = @(p,x) p(2).*x + p(1);
    p0 = [0.1, 1.0];
    fit_nonlin = fitnlm(log(Q_rec),log(-dQdt_rec),linFcn,p0,'Weight',weights);
    a = exp(fit_nonlin.Coefficients.Estimate(1));
    b = fit_nonlin.Coefficients.Estimate(2);
    
elseif strcmp(fitting_type, 'nonlinear') % nonlinear regression in lin space
    
    powerFcn = @(p,x) p(1).*x.^p(2);
    p0 = [0.1, 1.0];
    fit_nonlin = fitnlm(Q_rec,-dQdt_rec,powerFcn,p0,'Weight',weights);
    a = fit_nonlin.Coefficients.Estimate(1);
    b = fit_nonlin.Coefficients.Estimate(2);   
    
elseif strcmp(fitting_type, 'slope1') % linear regression in log log space with fixed slope 1
    
    linFcn = @(p,x) 1.*x + p(1);
    p0 = [0.1];
    fit_nonlin = fitnlm(log(Q_rec),log(-dQdt_rec),linFcn,p0,'Weight',weights);
    a = exp(fit_nonlin.Coefficients.Estimate(1));
    b = 1;

elseif strcmp(fitting_type, 'slope2') % linear regression in log log space with fixed slope 2

    linFcn = @(p,x) 2.*x + p(1);
    p0 = [0.1];
    fit_nonlin = fitnlm(log(Q_rec),log(-dQdt_rec),linFcn,p0,'Weight',weights);
    a = exp(fit_nonlin.Coefficients.Estimate(1));
    b = 2;
    
else
    error('Invalid fitting type.')
end

end


