function [FDC_slope, error_flag, error_str, fig_handles] = sig_FDC_slope(Q, t, varargin)
%sig_FDC_slope calculates slope of flow duration curve (FDC).
%   Calculates slope of FDC between two normalised streamflow percentiles
%   (e.g. 33rd and 66th, see McMillan et al., 2017). Slope can be fitted in
%   linear or in log space. Note that percentiles are defined as exceedance
%   probabilities, i.e. the low-section corresponds to high values, e.g.
%   [0.8 1.0].
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   slope_range: range of FDC [perc_lower perc_upper] in which slope should
%       be calculated, default is 33th to 66th percentile [0.33 0.66]
%   fitLogSpace: fit slope in log space or in linear space, default = true
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   FDC_slope: slope of flow duration curve [-]
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
%   FDC_slope = sig_FDC_slope(Q,t);
%   FDC_slope = sig_FDC_slope(Q,t,'slope_range',[0.33 0.66],'plot_results',true,'fitLogSpace',false);
%
%   References
%	McMillan, H., Westerberg, I. and Branger, F., 2017. Five guidelines for
%   selecting hydrological signatures. Hydrological Processes, 31(26),
%   pp.4757-4761.
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
addParameter(ip, 'slope_range', [0.33 0.66], validationFcn) % range of FDC
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
    FDC_slope = NaN;
    return
end

if slope_range(1)>slope_range(2) % check order of percentiles
    tmp = slope_range(1);
    slope_range(1) = slope_range(2);
    slope_range(2) = tmp;
end

% calculate signature
% get ranks as a proxy for exceedance probabilities
Q_sorted = sort(Q);
% Q_ranked = tiedrank(Q_sorted);
Q_ranked = [1:length(Q)]'; % give unique (random) rank to every measurement
FDC = 1 - Q_ranked./length(Q_ranked); % flow duration curve with unique ranks
Q_median = median(Q,'omitnan'); %

% slope of FDC between upper and lower percentile
indices = 1:length(FDC);

% due to the way the FDC is calculated, the maximum exceedance probability will always be <1
if slope_range(2) == 1
    bound_upp = 1; 
else
    bound_upp = max(indices(FDC >= slope_range(2)));
end
bound_low = max(indices(FDC >= slope_range(1)));

% fit slope, either in linear or in log space
if fitLogSpace
    FDC_slope = (log(Q_sorted(bound_upp)/Q_median) - log(Q_sorted(bound_low)/Q_median))./...
        (FDC(bound_upp) - FDC(bound_low));
else
    FDC_slope = ((Q_sorted(bound_upp)/Q_median) - (Q_sorted(bound_low)/Q_median))./...
        (FDC(bound_upp) - FDC(bound_low));
end

% in case flow is very intermittent (e.g. 66th percentile is 0)
if ~isfinite(FDC_slope)
    FDC_slope = NaN;
    error_flag = 3;
    error_str = ['Error: FDC slope could not be calculated, probably because flow is intermittent. ', error_str];
    return
end

% optional plotting
if plot_results
    fig = figure('pos',[100 100 350 300]); hold on
    if fitLogSpace
        plot(FDC,log(Q_sorted./Q_median),'linewidth',1.0)
        x = FDC(bound_low):0.001:FDC(bound_upp);
        c = log(Q_sorted(bound_low)./Q_median) - FDC_slope.*FDC(bound_low);
        y = FDC_slope.*x + c;
        ylabel('log(Q/median(Q)) [-]')
    else
        plot(FDC,(Q_sorted./Q_median),'linewidth',1.0)
        x = FDC(bound_low):0.001:FDC(bound_upp);
        c = (Q_sorted(bound_low)./Q_median) - FDC_slope.*FDC(bound_low);
        y = FDC_slope.*x + c;
        ylabel('Q/median(Q) [-]')
    end
    plot(x,y,'--','linewidth',2)
    legend('FDC','Fitted slope')
    xlabel('Exceedance probability [-]')
    fig_handles.FDC_slope = fig;
end

end

