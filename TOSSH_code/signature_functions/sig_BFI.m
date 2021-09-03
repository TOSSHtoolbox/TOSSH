function [BFI, error_flag, error_str, fig_handles] = sig_BFI(Q, t, varargin)
%sig_BFI calculates baseflow index (BFI).
%   Calculates BFI, that is the ratio between baseflow (volume) and
%   total streamflow (volume), with different baseflow separation methods 
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
%       filter_parameter assumed to be in days and converted accordingly
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   BFI: baseflow index [-]
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
%   BFI = sig_BFI(Q,t);
%   BFI_LH = sig_BFI(Q,t,'method','Lyne_Hollick','parameters',[0.925 3]);
%   BFI_UKIH = sig_BFI(Q,t,'method','UKIH','parameters',[5]);
%
%   References
%   Lyne, V. and Hollick, M., 1979. Stochastic time-variable
%   rainfall-runoff modelling. In Institute of Engineers Australia National
%   Conference (Vol. 1979, pp. 89-93). Barton, Australia: Institute of
%   Engineers Australia.
%   UK Institute of Hydrology (Great Britain), 1980. Low Flow Studies 
%   Reports. Institute of Hydrology.
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
    BFI = NaN;
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
            parameters = [filter_parameter, 3];
        elseif length(parameters) == 1
            parameters(1) = exp(log(parameters(1))/timestep_factor);
            parameters(2) = 3;
        elseif length(parameters) == 2
            parameters(1) = exp(log(parameters(1))/timestep_factor);        
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
        elseif length(parameters) == 1
            parameters(1) = parameters(1)*timestep_factor;
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

% calculate BFI
BFI = sum(Q_b,'omitnan')/sum(Q,'omitnan');

% check if 0<=BFI<=1
if BFI<0 || BFI>1   
    BFI = NaN; 
    error_flag = 1;
    error_str = ['Warning: Estimated BFI outside allowed range (0 to 1).', error_str];
end

% optional plotting
if plot_results
    fig = figure('Position',[100 100 700 250]); hold on
    plot(t,Q)
    plot(t,Q_b)
    xlabel('Date')
    ylabel('Flow [mm/timestep]')
    legend('Streamflow','Estimated baseflow')
    fig_handles.BFI = fig;
end

end
