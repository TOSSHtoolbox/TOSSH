function [PQ_strength, slope1, slope2, breakpoint, error_flag, error_str, fig_handles] ...
    = sig_PQ_Curve(Q, t, P, varargin)
%sig_PQ_Curve calculates signatures from cumulative P-Q curve.
%   Calculates cumulative differences between P and Q regime curves (i.e.
%   their averaged on each calendar day) and fits two straight lines 
%   separated by a breakpoint to it (broken stick fit). The slopes of these
%   lines, the breakpoint, as well the the ratio between the slopes
%   are returned as signatures (see Horner, 2020).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   OPTIONAL
%   interval: days of (water) year over which the two segments are fitted,
%       default = [15 183]
%   start_month: starting month, default = 10 (October)
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   PQ_strength: one minus ratio of slopes (0 if slopes are the same,
%       positive if slope 1 is larger and negative if slope 2 is larger)
%   slope1: slope of first segment ("dry slope")
%   slope2: slope of second segment ("wet slope")
%   breakpoint: date of breakpoint that defines the two segments
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
%   P = data.P;
%   [PQ_strength, slope1, slope2, breakpoint] = sig_PQ_Curve(Q,t,P);
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
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'P', @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1))

% optional input arguments
addParameter(ip, 'interval', [15 183], @(interval) isnumeric && size(interval)==2) % approx. first half year
addParameter(ip, 'start_month', 10, @(start_month) isnumeric(start_month) && numel(start_month)==1)
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, P, varargin{:})
interval = ip.Results.interval;
start_month = ip.Results.start_month;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P);
if error_flag == 2
    PQ_strength = NaN;
    slope1 = NaN;
    slope2 = NaN;
    breakpoint = NaN;
    return
end
timestep_days = days(timestep); % adjust for timestep

if any(interval) < 1 || any(interval > 365)
    error('Interval has to consist of values between 1 and 365.')
end

if any(start_month<1 | start_month>12) || any(floor(start_month)~=start_month)
    error('Month has to be a vector containing integers between 1 and 12.')
end

% calculate signature
% get average year
[Q_avg, t_avg] = util_AverageYear(Q,t,'start_water_year',start_month);
[P_avg, ~] = util_AverageYear(P,t,'start_water_year',start_month);

% calculate maximum difference between cumulative sums
P_cumsum = cumsum(P_avg,'omitnan')./timestep_days;
Q_cumsum = cumsum(Q_avg,'omitnan')./timestep_days;
PQ_diff = P_cumsum - Q_cumsum;
PQ_diff = PQ_diff(interval(1):interval(2))-PQ_diff(interval(1));

% fit two linear segments to get break point and slopes
% Horner (2020) forces first line to start at 0
xdata = [interval(1):interval(2)]'-interval(1);
ydata = PQ_diff;
dx = max(xdata) - min(xdata);
[breakpoint] = fminbnd(@(b2) util_FitBrokenStick(b2,xdata,ydata,true), ...
    min(xdata)+dx/100, max(xdata)-dx/100);
[err_b1,fittedlines,slopes] = util_FitBrokenStick(breakpoint,xdata,ydata,true);

% extract signatures
slope1 = slopes(1);
slope2 = slopes(2);
breakpoint = round(fittedlines(2,1)) + interval(1);
PQ_strength = 1 - slope2/slope1;

% optional plotting
if plot_results
    fig = figure('Position',[100 100 350 300]); hold on;
    p1=plot(xdata+ interval(1),ydata,'b','linewidth',1.5);
    p2=plot(fittedlines(1:2,1)+ interval(1),fittedlines(1:2,2),'r','linewidth',1.5);
    plot(fittedlines(2:3,1)+ interval(1),fittedlines(2:3,2),'r','linewidth',1.5)
    xlabel(strcat('Day of water year (starting in month',{' '},num2str(start_month),')'))
    ylabel('Cumulative P - Q [mm/timestep]')
    legend([p1,p2],{'Data','Fit'},'location','nw');
    fig_handles.PQ_curve = fig;
end

end

