function [Q_n_day_max, error_flag, error_str] = sig_Q_n_day_max(Q, t, n_day)
%Q_n_day_max calculates n day maxmimum of flow time series.
%   
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   n_day: window over which maximum should be calculated
%
%   OUTPUT
%   Q_n_day_max: n day maxmimum of flow [mm/n_day]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   Q_n_day_max = sig_Q_n_day_max(Q, t, 7);
%   Q_n_day_max = sig_Q_n_day_max(Q, t, [1:7]);
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
% n_day has to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'n_day', @(n_day) isnumeric(n_day) && (size(n_day,1)==1 || size(n_day,2)==1))

parse(ip, Q, t, n_day)

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    Q_n_day_max = NaN;
    return
end

if any(n_day<1 | n_day>length(t)) || any(floor(n_day)~=n_day)
    error('Month has to be a vector containing integers between 1 and length(timeseries).')
end

% calculate signature
Q_n_day_max = NaN(length(n_day),1);
for i=1:length(n_day)
    Q_n_day_max(i) = n_day(i)*max(movmean(Q,n_day(i))); % returns value in mm per chosen window n_day
end

end
