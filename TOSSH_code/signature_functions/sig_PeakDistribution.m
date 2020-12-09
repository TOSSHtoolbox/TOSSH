function [PeakDistribution, error_flag, error_str, fig_handles] = ...
    sig_PeakDistribution(Q, t, varargin)
%sig_PeakDistribution calculates peak distribution.
%   Calculates the the slope between the 10th and 50th of a flow duration 
%   curve constructed constructed by only considering hydrograph peaks 
%   (Euser et al., 2013). Slope can be fitted in linear or in log space.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   slope_range: range of peak FDC [perc_lower perc_upper] in which slope
%       should be calculated, default is 10th to 50th percentile [0.1 0.5]
%   fitLogSpace: fit slope in log space or in linear space, default = true
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   PeakDistribution: peak distribution [-]
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
%   PeakDistribution = sig_PeakDistribution(Q,t,'plot_results',true);
%   PeakDistribution = sig_PeakDistribution(Q,t,'plot_results',true,'fitLogSpace',false);
%
%   References
%   Euser, T., Winsemius, H.C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S.
%   and Savenije, H.H.G., 2013. A framework to assess the realism of model
%   structures using hydrological signatures. Hydrology and Earth System
%   Sciences, 17 (5), 2013.
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
validationFcn = @(x) isnumeric(x) && numel(x)==2 && (all(x >= 0)) && (all(x <= 1));
addParameter(ip, 'slope_range', [0.1 0.5], validationFcn) % range of FDC
addParameter(ip, 'fitLogSpace', true, @islogical) % fit slope in log space?
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, varargin{:})
slope_range = ip.Results.slope_range;
fitLogSpace = ip.Results.fitLogSpace;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    PeakDistribution = NaN;
    return
end

if slope_range(1)>slope_range(2) % check order of percentiles
    tmp = slope_range(1);
    slope_range(1) = slope_range(2);
    slope_range(2) = tmp;
end

% calculate signature
% get peaks
% Q_peak = findpeaks(Q);
Q_peak = Q(islocalmax(Q));


% get ranks as a proxy for exceedance probabilities
Q_peak_sorted = sort(Q_peak);
% Q_ranked = tiedrank(Q_sorted);
Q_peak_ranked = [1:length(Q_peak)]'; % give unique (random) rank to every measurement
FDC = 1 - Q_peak_ranked./length(Q_peak_ranked); % flow duration curve
Q_peak_median = median(Q_peak,'omitnan'); %

% slope of FDC between upper and lower percentile
indices = 1:length(FDC);
bound_upp = max(indices(FDC >= slope_range(2)));
bound_low = max(indices(FDC >= slope_range(1)));

% fit slope, either in linear or in log space
if fitLogSpace
    PeakDistribution = (log(Q_peak_sorted(bound_upp)/Q_peak_median) - ...
        log(Q_peak_sorted(bound_low)/Q_peak_median))./...
        (FDC(bound_upp) - FDC(bound_low));
else
    PeakDistribution = ((Q_peak_sorted(bound_upp)/Q_peak_median) - ...
        (Q_peak_sorted(bound_low)/Q_peak_median))./...
        (FDC(bound_upp) - FDC(bound_low));
end

% in case flow is very intermittent (e.g. 66th percentile is 0)
if ~isfinite(PeakDistribution)
    PeakDistribution = NaN;
    error_flag = 3;
    error_str = ['Error: Peak distibution could not be calculated. ', error_str];
    return
elseif isempty(PeakDistribution)     
    PeakDistribution = NaN;
    error_flag = 3;
    error_str = ['Error: No peaks detected. ', error_str];
    return    
end

% optional plotting
if plot_results
    fig = figure('pos',[100 100 350 300]); hold on
    hold on
    if fitLogSpace
        plot(FDC,log(Q_peak_sorted./Q_peak_median),'linewidth',1.0)
        x = FDC(bound_low):0.1:FDC(bound_upp);
        c = log(Q_peak_sorted(bound_low)./Q_peak_median) - PeakDistribution.*FDC(bound_low);
        y = PeakDistribution.*x + c;
        ylabel('log(Q/median(Q)) [-]')
    else
        plot(FDC,(Q_peak_sorted./Q_peak_median),'linewidth',1.0)
        x = FDC(bound_low):0.1:FDC(bound_upp);
        c = (Q_peak_sorted(bound_low)./Q_peak_median) - PeakDistribution.*FDC(bound_low);
        y = PeakDistribution.*x + c;
        ylabel('Q/median(Q) [-]')
    end
    plot(x,y,'--','linewidth',2)
    xlabel('Exceedance probability [-]')
    fig_handles.PeakDistribution = fig;
end

end

