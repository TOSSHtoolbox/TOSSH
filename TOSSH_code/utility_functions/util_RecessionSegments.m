function [flow_section, error_flag, error_str, fig_handles] = ...
    util_RecessionSegments(Q, t, varargin)
%util_RecessionSegments identifies all individual recession segments.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   recession_length: min. length of recessions (days), default = 5
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default = 0
%       (note that large eps values can lead to problematic recession
%       selection)
%   start_of_recession: define start of recession when baseflow filter
%       rejoins the curve "baseflow" or after peak "peak"
%   filter_par: smoothing parameter of Lyne-Hollick filter to determine
%      start of recession (higher = later recession start), default = 0.925
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   flow_section: n-by-2 array where n is the number of recession segments;
%       columns are the indices into the flow array of the start and end of
%       the recession segments
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
%   flow_section = util_RecessionSegments(Q,t);
%   flow_section = util_RecessionSegments(Q,t,'recession_length',5,...
%       'plot_results',true,'start_of_recession','peak','n_start',1);
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
addParameter(ip, 'recession_length', 5, @isnumeric) % length of decreasing
% flow in days to be declared a recession
addParameter(ip, 'n_start', 1, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'peak', @ischar) % defines start of a recession
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of
% Lyne-Hollick filter to determine start of recession (higher = later recession start)
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
plot_results = ip.Results.plot_results;
filter_par = ip.Results.filter_par;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
n_start = ip.Results.n_start;

% create empty figure handle
fig_handles = [];

% data checks

% default setting reads as good data
error_flag = 0;
error_str = '';

if eps > median(Q,'omitnan')/100
    error_flag = 1;
    error_str = ['Warning: eps set to a value larger than 1 percent of median(Q). High eps values can lead to problematic recession selection. ', error_str];
end

iszero = (Q==0);
Q(iszero) = NaN; % do not use days with zero flow

% identify all individual recession segments with length > recession_length days
% how many decreasing timesteps depends on length of timestep
len_decrease = recession_length/days(t(2)-t(1));
% find timesteps with decreasing flow
decreasing_flow = Q(2:end)<(Q(1:end-1)+eps);
% start on a non-decreasing point
start_point = find(decreasing_flow==0,1);
decreasing_flow = decreasing_flow(start_point:end);
% find start and end of decreasing sections
flow_change = find(decreasing_flow(1:end-1) ~= decreasing_flow(2:end));
% reshape into x by 2 array (columns = start, end of decrease)
flow_change = flow_change(1:(2*floor(size(flow_change,1)./2)));
flow_change = reshape(flow_change,2,[]).';
% find sections
flow_section = flow_change((flow_change(:,2)-flow_change(:,1))>=len_decrease+n_start,:);
flow_section = flow_section+start_point;
flow_section(:,1) = flow_section(:,1)+n_start; % move start point n days

% constantQ = false(length(flow_section),1);
% for i = 1:size(flow_section,1)
%     if length(unique(Q(flow_section(i,1):flow_section(i,2)))) == 1
%         constantQ(i) = true;
%     end
% end
% flow_section(constantQ,:) = [];

% % if eps > 0 it can happen that recessions have consecutive timesteps with
% % the same flow or 0 flow - these are removed
% if eps == 0
% else
%     rmv = false(length(flow_section),1);
%     for i = 1:size(flow_section,1)
%         Q_tmp = Q(flow_section(i,1):flow_section(i,2));
%         if any(diff(Q_tmp) == 0) || any(Q_tmp == 0) % diff diff
%             rmv(i) = true;
%         end
%     end
%     flow_section(rmv,:) = [];
% end

% remove recession segments that are just flat lines
% rmv = false(length(flow_section),1);
% for i = 1:size(flow_section,1)
%     Q_tmp = Q(flow_section(i,1):flow_section(i,2));
%     if max(Q_tmp)-min(Q_tmp) < 10^-8
%         rmv(i) = true;
%     end
% end
% flow_section(rmv,:) = [];

% % remove recession segments that contain flat lines of certain length
% rmv = false(length(flow_section),1);
% for i = 1:size(flow_section,1)
%     Q_tmp = Q(flow_section(i,1):flow_section(i,2));
%     k = 5;
%     if length(Q_tmp)>k
%         for j = 1:length(Q_tmp)-k
%             Q_tmp_sub = Q_tmp(j:j+k);
%             if max(Q_tmp_sub)-min(Q_tmp_sub) < 10^-6
%                 rmv(i) = true;
%             end
%         end
%     else
%         if max(Q_tmp)-min(Q_tmp) < 10^-6
%             rmv(i) = true;
%         end
%     end
% end
% flow_section(rmv,:) = [];

Q(iszero) = 0;

% optional plotting
if plot_results
    fig = figure('Position',[100 100 700 250]); hold on;
    title('Selected recession segments')
    h0 = plot(t,Q,'k','linewidth',1.5);
    for i = 1:size(flow_section,1)
        h1 = plot(t(flow_section(i,1):flow_section(i,2)),...
            Q(flow_section(i,1):flow_section(i,2)),'r-','linewidth',1.5);
        % plot(t(flow_section(i,1):flow_section(i,2)),Qm(flow_section(i,1):flow_section(i,2)),'r');
        % plot(t(flow_section(i,1):flow_section(i,2)),dQdt(flow_section(i,1):flow_section(i,2)),'b');
    end
    legend('Streamflow','Recessions')
    xlabel('Date')
    ylabel('Flow [mm/timestep]')
    fig_handles.RecessionSegments = fig;
end

switch start_of_recession
    
    case 'peak'        
        if size(flow_section,1)==0
            error_flag = 3;
            error_str = ['Error: No long enough recession periods, consider setting eps parameter > 0. ', error_str];
            return
        elseif size(flow_section,1) < 10
            error_flag = 1;
            error_str = ['Warning: Fewer than 10 recession segments extracted, results might not be robust. ', error_str];
        end
        
    case 'baseflow'
        % beginning of recession = point where baseflow filter rejoins the curve
        % (alternative needs basin area)
        % use Lyne-Hollick Filter with default a, nr_passes = 1
        % (doesn't work with nr_passes > 1 as then baseflow is never equal to flow)
        Q_b = util_LyneHollickFilter(Q,'filter_parameter',filter_par,'nr_passes',1);
        
        % find sections identified as baseflow in complete series
        isbaseflow = Q_b == Q;
        for i = 1:size(flow_section,1)
            % find which points in the decreasing section count as baseflow only
            isb_section = isbaseflow(flow_section(i,1):flow_section(i,2));
            if all(isb_section==0)
                flow_section(i,:) = NaN;
            else
                % find the first such point
                isb_start = find(isb_section==1,1,'first');
                % limit section to only those points
                flow_section(i,1)=flow_section(i,1)+isb_start-1;
            end
        end
        % if flow_section less than 4 points, remove
        flow_section = flow_section(~isnan(flow_section(:,1)),:);
        flow_section = flow_section(flow_section(:,2)>=flow_section(:,1)+3,:);
        if size(flow_section,1)==0
            error_flag = 3;
            error_str = ['Error: No long enough baseflow recession periods, consider increasing filter_par parameter. ', error_str];
            return
        elseif size(flow_section,1) < 10
            error_flag = 1;
            error_str = ['Warning: Fewer than 10 recession segments extracted, results might not be robust. ', error_str];
        end
        
        % add the baseflow recession sections to the plot
        if plot_results
            % add baseflow to plot
            h3 = plot(t,Q_b, 'b-','linewidth',1.5);
            for i = 1:size(flow_section,1)
                h2=plot(t(flow_section(i,1):flow_section(i,2)),...
                    Q(flow_section(i,1):flow_section(i,2)),'g-','linewidth',1.5);
            end
            legend([h0 h3 h1 h2],...
                {'Streamflow','Baseflow','Complete recessions', 'Baseflow recessions'})
            fig_handles.RecessionSegments = fig;
        end
        
    otherwise
        error('Incorrect option for start of recession.')
end

end