function [RecessionK_part, error_flag, error_str, fig_handles] = ...
    sig_RecessionParts(Q, t, part, varargin)
%sig_RecessionParts calculates recession constant for certain recession parts.
%   Fits exponential function (Q(t) = Q(0)*exp(-k*t)) to certain recession
%   parts (e.g. early or late recessions; see Horner, 2020). Recessions
%   are extracted from smoothed streamflow time series (moving mean),
%   recession maxima need to be larger than a threshold (e.g. median flow),
%   and recessions need to have a minimum length (default 5 days).
%   Early recessions are from day 1 to day 5 of a recession. Late
%   recessions are from day 15 to day 30 of a recession.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   part: recession part, e.g early, late, or custom
%   OPTIONAL
%   recession_part: recession part for which recession constant is
%       calculated [recession_start min_recession_length
%       max_recession_length]
%   Q_threshold: threshold above which recession maxima have to be,
%       default = median(Q)
%   moving_window: moving window used to smooth time series for recession
%       selection, default = 3
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   RecessionK_part: median recession constant for certain recession parts 
%       (e.g. early recessions)
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
%   k_median_early = sig_RecessionParts(Q,t,'early','plot_results',true);
%   k_median_late = sig_RecessionParts(Q,t,'late','plot_results',true);
%   k_median_early_custom = sig_RecessionParts(Q,t,'custom','recession_part',[1 5 5]);
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
% part has to be char and only one word
addRequired(ip, 'part', @(part) ischar(part) && size(part,1)==1)

% optional input arguments
addParameter(ip, 'recession_part', [], @isnumeric) % flow threshold
addParameter(ip, 'Q_threshold', median(Q,'omitnan'), @isnumeric) % Q threshold below which recessions are discarded
addParameter(ip, 'moving_window', 3, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, part, varargin{:})
recession_part = ip.Results.recession_part;
Q_threshold = ip.Results.Q_threshold;
moving_window = ip.Results.moving_window;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    RecessionK_part = NaN;
    return
end

% specify recession part
switch part
    case 'early'
        % early recessions are from day 1 to day 5 of a recession
        recession_start = 1;
        min_recession_length = 5;
        max_recession_length = 5;
        
    case 'late'
        % late recessions are from day 15 to day 30 of a recession and have
        % to be at least 5 days long
        recession_start = 15;
        min_recession_length = 5;
        max_recession_length = 15;
        
    case 'custom'
        if isempty(recession_part) || numel(recession_part) ~= 3
            error('No/wrong custom recession parts specified.')
        end
        recession_start = recession_part(1);
        min_recession_length = recession_part(2);
        max_recession_length = recession_part(3);
        
    otherwise
        error('Incorrect recession part specified.')
end

% calculate signature

% get recession segments
Q_tmp = movmean(Q,moving_window,'omitnan'); % smooth time series before recession extraction
Q_tmp(isnan(Q)) = NaN;
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[flow_section, error_flag, error_str, fig_handles] = ...
    util_RecessionSegments(Q_tmp, t, 'plot_results', plot_results, ...
    'recession_length', min_recession_length, 'n_start', recession_start - 1); % discard all values before recession start
if error_flag == 3
    RecessionK_part = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

% remove recessions that have too low streamflow values
below_threshold = false(size(flow_section,1),1);
for i = 1:size(flow_section,1)
    rec = [flow_section(i,1):flow_section(i,2)]'; % get recession
    if max(Q(rec)) < Q_threshold
        below_threshold(i) = true;
    end
end
flow_section(below_threshold,:) = [];

% truncate recessions if they are longer than maximum recession lengths
for i = 1:size(flow_section,1)
    if flow_section(i,1) + max_recession_length - 1 < flow_section(i,2)
        flow_section(i,2) = flow_section(i,1) + max_recession_length - 1;
    end
end

if isempty(flow_section)
    RecessionK_part = NaN;
    error_flag = 3;
    error_str = ['Error: No recession periods that match the recession selection criteria. ', error_str];
    return
end

% fit exponential to recessions
k = NaN(size(flow_section,1),1);
for i = 1:size(flow_section,1)
    rec = [flow_section(i,1):flow_section(i,2)]'; % get recession
    % default: linear fit in log log space for fast calculation
    k(i) = util_FitExponential(Q(rec),t(rec),'semilog'); 
end

% get median recession constant
RecessionK_part = median(k,'omitnan');

% optional plotting
if plot_results
    fig = figure('Position',[100 100 700 250]); hold on
    p1 = plot(t,Q,'k','linewidth',1.5);
    for i = 1:size(flow_section,1)
        rec = [flow_section(i,1):flow_section(i,2)]'; % get recession
        t_rec = [0:length(rec)-1]';
        Q_est = Q(rec(1))*exp(t_rec.*-k(i));
        p2 = plot(t(rec),Q_est,'g','linewidth',1.5);
    end
    xlabel('Date')
    ylabel('Flow [mm/timestep]')
    title('Exponential recessions fitted to certain recession parts')
    legend([p1 p2],{'Streamflow','Fitted recessions'},'Location','best')
    % set(gca,'yscale','log');
    fig_handles.RecessionParts = fig;
end

end

