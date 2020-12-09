function [AC, error_flag, error_str] = sig_Autocorrelation(Q, t, varargin)
%sig_Autocorrelation caculates lag-x autocorrelation of flow.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   lag: time lag at which autocorrelation should be calculated (default =
%       1 timestep), can also be a vector
%
%   OUTPUT
%   AC: lag-x autocorrelation [-], can also be a vector
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   AC = sig_Autocorrelation(Q,t);
%   AC = sig_Autocorrelation(Q,t,'lag',1);
%   AC = sig_Autocorrelation(Q,t,'lag',[0:length(Q)-2]);
%
%   References
%   https://en.wikipedia.org/wiki/Autocorrelation
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))

% optional input arguments
addParameter(ip, 'lag', 1, @isnumeric) % lag for autocorrelation

parse(ip, Q, t, varargin{:})
lag = ip.Results.lag;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    AC = NaN;
    return
end

if any(lag>length(Q)-2) || any(lag<0)
    error('lag must be between 0 and must not exceed the number of observations minus two.')
end

% calculate signature (different methods use different toolboxes)
AC = nanxcov(Q,Q,max(lag),'coeff');
AC = AC(max(lag)+1+lag);

% AC = autocorr(Q,max(lag)); % from Econometrics toolbox
% AC = AC(lag+1);

end
