function [x_Q_frequency, error_flag, error_str] = sig_x_Q_frequency(Q, t, type, varargin)
%sig_x_Q_frequency calculates various kinds of flow frequencies.
%   Calculates various kinds of flow frequencies, e.g. no flow frequency, 
%   high flow frequency. Typical metrics can be chosen from standard list 
%   (no, high, low), or created manually (e.g. 0.5*median(Q)).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   type: type of flow frequency (no, high, low, custom_high, custom_low)
%   OPTIONAL
%   threshold: flow threshold above (below) flow frequency is calculated 
%       (e.g. 9*median(Q) for high) [mm/timestep]
%
%   OUTPUT
%   x_Q_freq: x flow frequency [-]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   x_Q_freq = sig_x_Q_freq(Q,t,'no');
%   x_Q_freq = sig_x_Q_freq(Q,t,'high');
%   x_Q_freq = sig_x_Q_freq(Q,t,'low');
%   x_Q_freq = sig_x_Q_freq(Q,t,'custom_low','threshold',0.2*mean(Q,'omitnan'));
%
%   References
%   Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and  
%   Clark, M.P., 2018. A ranking of hydrological signatures based on their 
%   predictability in space. Water Resources Research, 54(11), pp.8792-8812.
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
% type has to be char and only one word
addRequired(ip, 'type', @(type) ischar(type) && size(type,1)==1) 

% optional input arguments
addParameter(ip, 'threshold', [], @isnumeric) % flow threshold

parse(ip, Q, t, type, varargin{:})
threshold = ip.Results.threshold;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    x_Q_frequency = NaN;
    return
end

% calculate signature
switch type
    
    case 'no'
        x_Q_num = length(Q(Q==0));
        x_Q_frequency = x_Q_num/length(Q);
        
    case 'high'        
        Q_high = 9*median(Q,'omitnan');
        x_Q_num = length(Q(Q>Q_high));
        x_Q_frequency = x_Q_num/length(Q);
        
    case 'low'
        Q_low = 0.2*mean(Q,'omitnan');
        x_Q_num = length(Q(Q<Q_low));
        x_Q_frequency = x_Q_num/length(Q);
        
    case 'custom_high'
        if isempty(threshold) || numel(threshold) > 1
            error('No/wrong custom threshold specified.')
        end
        x_Q_num = length(Q(Q>threshold));
        x_Q_frequency = x_Q_num/length(Q);
        
    case 'custom_low'
        if isempty(threshold) || numel(threshold) > 1
            error('No/wrong custom threshold specified.')
        end
        x_Q_num = length(Q(Q<threshold));
        x_Q_frequency = x_Q_num/length(Q);
        
    otherwise
        error('Incorrect flow frequency type specified.')
end

end
