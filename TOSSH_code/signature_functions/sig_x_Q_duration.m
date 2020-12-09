function [x_Q_duration, error_flag, error_str] = sig_x_Q_duration(Q, t, type, varargin)
%sig_x_Q_duration calculates various kinds of flow durations.
%   Calculates various kinds of flow durations, e.g. no flow duration, 
%   high flow duration. Typical metrics can be chosen from a standard list
%   (no, high, low), or created manually (e.g. 0.5*median(Q)).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   type: type of flow duration (no, high, low, custom_high, custom_low)
%   OPTIONAL
%   threshold: flow threshold above (below) flow duration is calculated 
%       (e.g. 9*median(Q) for high) [mm/timestep]
%
%   OUTPUT
%   x_Q_dur: x flow duration [timestep]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   x_Q_dur = sig_x_Q_dur(Q,t,'no');
%   x_Q_dur = sig_x_Q_dur(Q,t,'high');
%   x_Q_dur = sig_x_Q_dur(Q,t,'low');
%   x_Q_dur = sig_x_Q_dur(Q,t,'custom_high','threshold',9*median(Q,'omitnan'));
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
    x_Q_duration = NaN;
    return
end

% calculate signature
switch type
    
    case 'no'
        no_Q = Q==0;         
        % find consecutive timesteps with no flows
        start1 = strfind([0,no_Q'],[0 1]);
        end1 = strfind([no_Q',0],[1 0]);
        interval_lengths = end1 - start1 + 1;
        
    case 'high'        
        Q_high = 9*median(Q,'omitnan');
        high_Q = Q>Q_high;
        % find consecutive timesteps with high flows
        start1 = strfind([0,high_Q'],[0 1]);
        end1 = strfind([high_Q',0],[1 0]);
        interval_lengths = end1 - start1 + 1;
        
    case 'low'
        Q_low = 0.2*mean(Q,'omitnan');
        low_Q = Q<Q_low; 
        % find consecutive timesteps with low flows
        start1 = strfind([0,low_Q'],[0 1]);
        end1 = strfind([low_Q',0],[1 0]);
        interval_lengths = end1 - start1 + 1;
        
    case 'custom_high'
        if isempty(threshold) || numel(threshold) > 1
            error('No/wrong custom threshold specified.')
        end
        custom_Q = Q>threshold; 
        % find consecutive timesteps above x flows
        start1 = strfind([0,custom_Q'],[0 1]);
        end1 = strfind([custom_Q',0],[1 0]);
        interval_lengths = end1 - start1 + 1;
        
    case 'custom_low'
        if isempty(threshold) || numel(threshold) > 1
            error('No/wrong custom threshold specified.')
        end
        custom_Q = Q<threshold; 
        % find consecutive timesteps below x flows
        start1 = strfind([0,custom_Q'],[0 1]);
        end1 = strfind([custom_Q',0],[1 0]);
        interval_lengths = end1 - start1 + 1;
        
    otherwise
        error('Incorrect flow duration type specified.')
end

if isempty(interval_lengths)
    x_Q_duration = 0;
else
    x_Q_duration = mean(interval_lengths);
end

end
