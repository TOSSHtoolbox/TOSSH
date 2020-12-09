function [SnowDayRatio, error_flag, error_str] = sig_SnowDayRatio(Q, t, P, T, varargin)
%sig_SnowDayRatio calculates snow day ratio.
%   The snow day ratio is defined as the number of days that experience
%   precipitation when the average daily air temperature is below 2 degC, 
%   divided by the total number of days per year with precipitation (see 
%   e.g. Sawicz et al., 2011).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   T: temperature [degC]
%   OPTIONAL
%   T_threshold: temperature threshold, default = 2 degC
%   
%   OUTPUT
%   SnowDayRatio: snow day ratio [-]
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
%   T = data.T;
%   SnowDayRatio = sig_SnowDayRatio(Q,t,P,T);
%
%   References
%   Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G.,
%   2011. Catchment classification: empirical analysis of hydrologic
%   similarity based on catchment function in the eastern USA. Hydrology
%   and Earth System Sciences, 15(9), pp.2895-2911.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 4
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
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'T', @(T) isnumeric(T) && (size(T,1)==1 || size(T,2)==1)) 

% optional input arguments
addParameter(ip, 'T_threshold', 2, @isnumeric) % temperature threshold

parse(ip, Q, t, P, T, varargin{:})
T_threshold = ip.Results.T_threshold;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P, 'T', T);
if error_flag == 2
    SnowDayRatio = NaN;
    return
end

% calculate signature
if numel(T_threshold)~=1
    error('Temperature threshold has to be a single number.')
elseif T_threshold<-10 || T_threshold > 10
    warning('Threshold should be somewhere around 0°C.')
end

NP = sum(P>0);
NS = sum(P>0 & T<T_threshold);
SnowDayRatio = NS/NP;

end
