function [FlashinessIndex, error_flag, error_str] = sig_FlashinessIndex(Q, t, varargin)
%sig_FlashinessIndex calculates Richards-Baker flashiness index.
%   Flashiness index defined as sum of absolute differences between 
%   consecutive (daily) flows following Baker et al. (2004).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%
%   OUTPUT
%   FlashinessIndex: Richards-Baker flashiness index [-]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;  
%   FlashinessIndex = sig_FlashinessIndex(Q,t);
%
%   References
%   Baker, D.B., Richards, R.P., Loftus, T.T. and Kramer, J.W., 2004. A new 
%   flashiness index: Characteristics and applications to midwestern rivers
%   and streams 1. JAWRA Journal of the American Water Resources 
%   Association, 40(2), pp.503-522.
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

parse(ip, Q, t, varargin{:})

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    FlashinessIndex = NaN;
    return
end

% calculate signature
FlashinessIndex = sum(abs(diff(Q)),'omitnan')/sum(Q,'omitnan');

end

