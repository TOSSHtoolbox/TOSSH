function [BaseflowMagnitude, error_flag, error_str, fig_handles] = ...
    sig_BaseflowMagnitude(Q, t, varargin)
%sig_BaseflowMagnitude calculates baseflow regime magnitude.
%   Calculates the difference between the minimum and the maximum of the
%   baseflow regime, defined as the average baseflow on each calendar day 
%   (see Horner, 2020). Different baseflow separation methods can be used 
%   (Lyne and Hollick, 1979; UK Institute of Hydrology, 1980).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   method: which baseflow separation method should be employed 
%       ('Lyne_Hollick','UKIH'), default = 'UKIH'
%   parameters: specify filter parameters ([filter_parameter nr_passes] for 
%       Lyne Hollick, default = [0.925 3]; [n_days] for UKIH, default = 5)
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   BaseflowMagnitude: baseflow regime magnitude [mm]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%   fig_handles: figure handles to manipulate figures (empty if plotting is
%       not requested)
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   BaseflowMagnitude = sig_BaseflowMagnitude(Q,t);
%
%   References
%   Horner, I., 2020. Design and evaluation of hydrological signatures for 
%   the diagnostic and improvement of a process-based distributed 
%   hydrological model (Doctoral dissertation, Universite Grenoble Alpes).
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
addParameter(ip, 'method', 'UKIH', @ischar) % which method? Default: UKIH
addParameter(ip, 'parameters', [], @isnumeric) % which parameter values?
% addParameter(ip, 'threshold_type', [], @ischar) % how to threshold Lyne-Hollick filter?
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, varargin{:})
method = ip.Results.method;
parameters = ip.Results.parameters;
% specify when to threshold (default: after each pass)
threshold_type = 'pass'; % ip.Results.threshold_type;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    BaseflowMagnitude = NaN;
    return
end
timestep_factor = 1/days(timestep); % adjust for timestep

% calculate signature

% pad time series to compensate for warm up effect (Ladson et al., 2013)
if length(Q)>60
    Q_padded = [Q(30:-1:1); Q; Q(end-29:end)];
else
    Q_padded = Q;
    error_flag = 1;
    error_str = ['Warning: Very short time series. Baseflow separation might be unreliable. ', error_str];
end

% obtain baseflow
switch method
    
    case 'Lyne_Hollick'
        if isempty(parameters)
            filter_parameter = exp(log(0.925)/timestep_factor);
            parameters = [filter_parameter, 1];
        elseif length(parameters) == 1
            parameters(2) = 1;
        elseif length(parameters) > 2
            error('Too many filter parameters.')
        end
        
        if isempty(threshold_type)
            Q_b = util_LyneHollickFilter(Q_padded, ...
                'filter_parameter', parameters(1), 'nr_passes', parameters(2));
        else
            Q_b = util_LyneHollickFilter(Q_padded, ...
                'filter_parameter', parameters(1), 'nr_passes', parameters(2),...
                'threshold_type',threshold_type);
        end
        
    case 'UKIH'
        if isempty(parameters)
            parameters = 5*timestep_factor;
        elseif length(parameters) > 1
            error('Too many filter parameters.')
        end
        Q_b = util_UKIH_Method(Q_padded, 'n_days', parameters(1));
        
    otherwise
        error('Please choose one of the available baseflow separation methods (UKIH or Lyne_Hollick).')
end

% remove padding
if length(Q)>60
    Q_b(1:30) = [];
    Q_b(end-29:end) = [];
else
end

% calculate baseflow regime
[B_regime,t_regime,fig_handles] = util_AverageYear(Q_b,t,'plot_results',plot_results);

% calculate baseflow regime magnitude
BaseflowMagnitude = (max(B_regime) - min(B_regime))*timestep_factor;

end
