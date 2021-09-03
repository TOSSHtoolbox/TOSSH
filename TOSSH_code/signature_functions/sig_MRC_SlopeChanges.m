function [MRC_num_segments, Segment_slopes, error_flag, error_str, fig_handles] = ...
    sig_MRC_SlopeChanges(Q, t, varargin)
%sig_MRC_SlopeChanges calculates MRC and whether it contains significant changes in slope.
%   According to Estrany et al. (2010), Clark et al. (2009), and others,
%   the changes in slope represent different reservoirs contributing to the
%   runoff response.
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%	OPTIONAL
%   recession_length:  min. length of recessions [days], default = 15
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default = 0
%   start_of_recession: define start of recession when baseflow filter
%       rejoins the curve ("baseflow"), or after hydrograph peak ("peak")
%   filter_par: smoothing parameter of Lyne-Hollick filter to determine
%      start of recession (higher = later recession start), default = 0.925
%   seg_test: reduction in RMSE needed to recommend extra segment in MRC,
%       default = 0.75
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   MRC_num_segments: number of different segments in MRC
%   Segment_slopes: slopes of master recession curve segments in Q vs. 
%       relative time plot [1/timestep]
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
%   MRC_num_segments = sig_MRC_SlopeChanges(Q, t);
%
%   References
%   Estrany, J., Garcia, C. and Batalla, R.J., 2010. Hydrological response
%   of a small mediterranean agricultural catchment. Journal of Hydrology,
%   380(1-2), pp.180-190.
%   Clark, M.P., Rupp, D.E., Woods, R.A., Tromp-van Meerveld, H.J., Peters,
%   N.E. and Freer, J.E., 2009. Consistency between hydrological models and
%   field observations: linking processes at the hillslope scale to
%   hydrological responses at the watershed scale. Hydrological Processes,
%   23(2), pp.311-319.
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
addParameter(ip, 'ignoreNaN', 'y', @ischar) % ignore NaN values y/n?
addParameter(ip, 'recession_length', 15, @isnumeric) % length of decreasing
% flow in days to be declared a recession
addParameter(ip, 'n_start', 0, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'peak', @ischar) % defines start of a recession
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of
% Lyne-Hollick filter to determine start of recession (higher = later recession start)
addParameter(ip, 'seg_test', 0.75, @isnumeric) % what reduction in RMSE
% needed to recommend extra segment in MRC
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
n_start = ip.Results.n_start;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
filter_par = ip.Results.filter_par;
plot_results = ip.Results.plot_results;
seg_test = ip.Results.seg_test;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    MRC_num_segments = NaN;
    Segment_slopes = NaN;
    return
end

% calculate signature

% identify all individual recession segments
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[flow_section, error_flag, error_str, fig_handles] = util_RecessionSegments(Q, t, ...
    'recession_length', recession_length, 'eps', eps, ...
    'filter_par', filter_par, 'plot_results', plot_results, ...
    'start_of_recession', start_of_recession, 'n_start', n_start);
if error_flag == 3
    MRC_num_segments = NaN;
    Segment_slopes = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

% MRC constucted using the adapted matching strip method (Posavec et al., 2006)
[mrc] = util_MasterRecessionCurve(Q, flow_section, ...
    'fit_method', 'nonparametric_analytic', 'plot_results', false);

% k = slope of the linear regression between log-transformed discharge and recession length
mdl = [(mrc(:,1).^0) (mrc(:,1))]\log(mrc(:,2));
slope_b0 = mdl(2);
err_b0 = sqrt(sum((log(mrc(:,2))-(mdl(1) + mdl(2).*(sort(mrc(:,1))))).^2));

% optimise two-segment fit
xdata=mrc(:,1);
ydata=log(mrc(:,2));

dx = max(xdata) - min(xdata);
[breaks] = fminbnd(@(b2) util_FitBrokenStick(b2,xdata,ydata),...
    min(xdata) + dx/100, max(xdata) - dx/100);
[err_b1,fittedlines,slopes_b1] = util_FitBrokenStick(breaks,xdata,ydata);

% optimise three-segment fit
options = optimoptions('fmincon'); 
options.Display = 'off'; % to stop fmincon from displaying info
% get initialization points for breaks
bk1 = min(xdata) + dx/3;
bk2 = min(xdata) + dx*2/3;
% check that there are at least 10 points in the last section and revise breakpoints if not
if sum(xdata > bk2) < 10
    bk2 = xdata(end-10);
    bk1 = min(xdata) + (bk2 - min(xdata))/2;
end

[breaks2] = fmincon(@(b2) util_FitBrokenStick(b2,xdata,ydata),...
    [bk1, bk2],[],[],[],[],(min(xdata) + dx/100)*ones(2,1),...
    (max(xdata) - dx/100)*ones(2,1),[],options);
[err_b2,fittedlines2,slopes_b2] = util_FitBrokenStick(breaks2,xdata,ydata);

RMS_errors = [err_b0,err_b1,err_b2];

MRC_num_segments = 1;
Segment_slopes = slope_b0;
if RMS_errors(2)<seg_test*RMS_errors(1)
    MRC_num_segments=2;
    Segment_slopes = slopes_b1;
    if RMS_errors(3)<seg_test*RMS_errors(2)
        MRC_num_segments=3;
        Segment_slopes = slopes_b2;
    end
elseif RMS_errors(3)<(seg_test.^2)*RMS_errors(1)
    MRC_num_segments=3;
    Segment_slopes = slopes_b2;
end

Segment_slopes = -Segment_slopes; % convert slope to recession rate

% optional plotting
if plot_results
    % plot results to demonstrate fit
    fig = figure('Position',[100 100 350 300]);    
    semilogy(mrc(:,1),(mrc(:,2)),'k-','linewidth',2)
    hold on
    semilogy(sort(mrc(:,1)),exp(mdl(1) + mdl(2).*(sort(mrc(:,1)))),'b-')
    % plot broken stick fits
    semilogy(fittedlines(:,1),exp(fittedlines(:,2)),'g-')
    semilogy(fittedlines2(:,1),exp(fittedlines2(:,2)),'r-')    
    xlabel('Relative time [timestep]')
    ylabel('Flow [mm/timestep]')
    title('Fitted recession curve - log scale')
    legend('Master Recession Curve',['1-segment fit, RMSE = ',num2str(err_b0)],...
        ['2-segment fit, RMSE = ',num2str(err_b1)],...
        ['3-segment fit, RMSE = ',num2str(err_b2)],'location','best')
    fig_handles.MRC_slopeChanges = fig;
end

end
