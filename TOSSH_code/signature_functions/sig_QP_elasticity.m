function [QP_elasticity, error_flag, error_str] = sig_QP_elasticity(Q, t, P, varargin)
%sig_QP_elasticity calculates streamflow-precipitation elasticity.
%   Calculates streamflow-precipitation elasticity, the sensitivity of a
%   catchment's streamflow response to changes in precipitation at the
%   annual time scale. Start of water year can be defined.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   OPTIONAL
%   start_water_year: first month of water year, default = 10 (October)
%   method: choose method to calculate elasticity ('Sawicz','Sanka'),
%       default = 'Sanka'
%
%   OUTPUT
%   QP_elasticity: streamflow-precipitation elasticity [-]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   P = data.P;
%   QP_elasticity = sig_QP_elasticity(Q,t,P,'start_water_year',10);
%
%   References
%   Sankarasubramanian, A., Vogel, R.M. and Limbrunner, J.F., 2001. Climate
%   elasticity of streamflow in the United States. Water Resources
%   Research, 37(6), pp.1771-1781.
%	Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G.,
%   2011. Catchment classification: empirical analysis of hydrologic
%   similarity based on catchment function in the eastern USA. Hydrology
%   and Earth System Sciences, 15(9), pp.2895-2911.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 3
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'P', @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1))

% optional input arguments
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 1) && (x <= 12) && floor(x)==x;
addParameter(ip, 'start_water_year', 10, validationFcn) % When does the water year start? Default: 10
addParameter(ip, 'method', 'Sanka', @ischar) % which method? Default: Sanka

parse(ip, Q, t, P, varargin{:})
start_water_year = ip.Results.start_water_year;
method = ip.Results.method;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P);
if error_flag == 2
    QP_elasticity = NaN;
    return
end

% calculate signature
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
% aggregate time series to get annual sums
[Q_annual, ~, ~, error_flag, error_str] = util_AggregateTimeSeries(Q, t, start_water_year);
[P_annual, ~, ~, error_flag, error_str] = util_AggregateTimeSeries(P, t, start_water_year);
if error_flag == 0
    error_flag = error_flag_tmp;
    error_str = error_str_tmp;
end

% calculate elasticity
switch method
    case 'Sanka'
        dQ = Q_annual-mean(Q_annual,'omitnan');
        dP = P_annual-mean(P_annual,'omitnan');
    case 'Sawicz'
        dQ = diff(Q_annual);
        dP = diff(P_annual);
    otherwise
        error('Please choose one of the available baseflow separation methods (Sawicz or Sanka).')
end

QP_elasticity = median(...
    (dQ./dP)*(mean(P_annual,'omitnan')/mean(Q_annual,'omitnan')),'omitnan');

end
