function [Q_mean_monthly, error_flag, error_str] = sig_Q_mean_monthly(Q, t, month)
%sig_Q_mean_monthly calculates mean monthly flow of time series.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   month: month(s) in which mean flow should be calculated (can be vector)
%
%   OUTPUT
%   Q_mean_monthly: mean flow in specified month(s) [mm/timestep]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   Q_mean_monthly = sig_Q_mean_monthly(Q, t, 8); % August
%   Q_mean_monthly = sig_Q_mean_monthly(Q, t, [10:12 1:9]); % full year
%
%   References
%   https://en.wikipedia.org/wiki/Arithmetic_mean
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
% month has to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'month', @(month) isnumeric(month) && (size(month,1)==1 || size(month,2)==1)) 

parse(ip, Q, t, month)

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    Q_mean_monthly = NaN;
    return
end

if any(month<1 | month>12) || any(floor(month)~=month)
    error('Month has to be a vector containing integers between 1 and 12.')
end

% calculate signature
date_vec = datevec(t);
Q_mean_monthly = NaN(length(month),1);
for i = 1:length(month)
    Q_mean_monthly(i) = mean(Q(date_vec(:,2)==month(i)),'omitnan'); % always ignoring NaNs 
end

end
