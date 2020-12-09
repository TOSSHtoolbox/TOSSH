function [CoV, error_flag, error_str] = sig_Q_CoV(Q, t)
%sig_Q_CoV calculates coefficient of variation (CoV) of flow.
%   Calculates CoV, i.e. standard deviation normalised by mean.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%
%   OUTPUT
%   CoV: coefficient of variation [-]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   CoV = sig_Q_CoV(Q,t);
%
%   References
%   https://en.wikipedia.org/wiki/Coefficient_of_variation
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

parse(ip, Q, t)

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    CoV = NaN;
    return
end

% calculate signature
CoV = std(Q,'omitnan')/mean(Q,'omitnan');

end
