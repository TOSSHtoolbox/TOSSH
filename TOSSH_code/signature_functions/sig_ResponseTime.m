function [ResponseTime, error_flag, error_str] = sig_ResponseTime(Q, t, P, varargin)
%sig_ResponseTime calculates catchment response time.
%   Calculates the catchment response time given a rainfall and a
%   streamflow time series using the DCMA method from Giani et al. (2020).
%
%   Notes:
%   Works best for sub-daily data. 
%   Can be slow if time series are long (e.g. long hourly time series).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   OPTIONAL
%   max_window: Maximum window [timestep] tested. Set it sensibly according
%       to the resolution of your data (e.g. for hourly data,
%       max_window = 300 means that time of concentration can be maximum
%       300hours/2 = 150hours =~ 6days). Default is 15 days.
%
%   OUTPUT
%   ResponseTime: catchment response time [timestep]
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
%   ResponseTime = sig_ResponseTime(Q,t,P);
%
%   References
%   Giani et al., 2020. A Practical, Objective and Robust Technique to
%   Directly Estimate Catchment Response Time, submitted to WRR
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
addParameter(ip, 'max_window', [], @isnumeric) % maximum window [days]

parse(ip, Q, t, P, varargin{:})
max_window = ip.Results.max_window;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P);
if error_flag == 2
    ResponseTime = NaN;
    return
end

% adjust max_window to time step
if isempty(max_window)
    max_window = 15/days(timestep);
elseif max_window < 3 || max_window > length(Q)
    error('Window size cannot be smaller than 3 times the timestep or larger than length of time series.')
elseif  max_window*days(timestep) > 100
    error_flag = 1;
    error_str = ['Warning: Window size is very large. ', error_str];
end

% calculate signature
P = P';
Q = Q';
P_int = cumsum(P, 'omitnan'); % cumulating rainfall time series (Eq.1)
Q_int = cumsum(Q, 'omitnan'); % cumulating streamflow time series (Eq.2)
len = length(P); % length of the time series

% todo: intialise arrays?

for w=3:2:max_window
    P_mean((w-1)/2,:) = movmean(P_int, w); % moving average on the integrated rainfall time series (Eq.5)
    Q_mean((w-1)/2,:) = movmean(Q_int, w); % moving average on the integrated streamflow time series (Eq.6)
    
    flutt_P((w-1)/2,:) = P_int-P_mean((w-1)/2,:);
    F_P((w-1)/2) = (1/(len-w+1))*...
        sum((flutt_P((w-1)/2,w-0.5*(w-1):len-0.5*(w-1))).^2,'omitnan'); % squared rainfall fluctuations (Eq.3)
    
    flutt_Q((w-1)/2,:) = Q_int-Q_mean((w-1)/2,:);
    F_Q((w-1)/2) = (1/(len-w+1))*...
        sum((flutt_Q((w-1)/2,w-0.5*(w-1):len-0.5*(w-1))).^2,'omitnan'); % squared streamflow fluctuations (Eq.4)
    
    F_PQ((w-1)/2) = (1/(len-w+1))*...
        sum(flutt_P((w-1)/2,w-0.5*(w-1):len-0.5*(w-1)).*...
        flutt_Q((w-1)/2,w-0.5*(w-1):len-0.5*(w-1)),'omitnan'); % bivariate rainfall-streamflow fluctuations (Eq.7)
    rho((w-1)/2) = F_PQ((w-1)/2)/(sqrt(F_P((w-1)/2))*sqrt(F_Q((w-1)/2))); % DMCA-based correlation coefficent (Eq.8)
end

position_minimum = find(rho==min(rho,[],'omitnan'));
if isempty(position_minimum)
    ResponseTime = NaN;
    error_flag = 3;
    error_str = ['Error: Response time could not be calculated. ', error_str];
    return
else
    ResponseTime = position_minimum;
end

end