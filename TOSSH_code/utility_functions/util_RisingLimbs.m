function [flow_section, error_flag, error_str, fig_handles] = ...
    util_RisingLimbs(Q, t, varargin)
%util_RisingLimbs identifies all rising limbs.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   rising_limb_length: length of rising limbs [days], default = 1
%   eps: allowed decrease in flow during rising limb, default = 0
%   minimum_peak: minimum peak to be counted as rising limb (peak size is
%       defined as difference between end and start of rising limb)
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   flow_section: n-by-2 array where n is the number of rising limbs
%       columns are the indices into the flow array of the start and end of
%       the rising limbs
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
%   flow_section = util_RisingLimbs(Q, t);
%   flow_section = util_RisingLimbs(Q, t, 'rising_limb_length', 2, 'plot_results', true);
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
addParameter(ip, 'rising_limb_length', 1, @isnumeric) % length of increasing flow in days to be declared a rising limb
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during rising limb
addParameter(ip, 'minimum_peak', 0, @isnumeric) % minimum peak to be counted as rising limb
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
rising_limb_length = ip.Results.rising_limb_length;
eps = ip.Results.eps;
minimum_peak = ip.Results.minimum_peak;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% default setting reads as good data
error_flag = 0;
error_str = '';

% identify all individual rising limbs with length > rising_limb_length days
% how many increasing timesteps depends on length of timestep
len_increase = rising_limb_length/days(t(2)-t(1));
% hind timesteps with increasing flow
increasing_flow = Q(2:end)>(Q(1:end-1)-eps);
% start on a non-increasing point
start_point = find(increasing_flow==0,1);
increasing_flow = increasing_flow(start_point:end);
% find start and end of increasing sections
flow_change = find(increasing_flow(1:end-1) ~= increasing_flow(2:end));
% reshape into x by 2 array (columns = start, end of decrease)
flow_change = flow_change(1:(2*floor(size(flow_change,1)./2)));
flow_change = reshape(flow_change,2,[]).';
% find sections
flow_section = flow_change((flow_change(:,2)-flow_change(:,1))>=len_increase,:);
flow_section = flow_section+start_point;
flow_section(:,1) = flow_section(:,1); % move start point n days
% remove rising limbs which have a peak lower than minimum_peak
% flow_section((Q(flow_section(:,2)) < minimum_peak),:) = [];
% remove rising limbs which have a peak lower than minimum_peak
flow_section((Q(flow_section(:,2)) - Q(flow_section(:,1))) < minimum_peak,:) = [];

if numel(flow_section)==0
    error_flag = 3;
    error_str = ['Error: No long enough rising limbs, consider setting eps parameter > 0. ', error_str];
end

% optional plotting
if plot_results
    fig = figure('Position',[100 100 700 250]); hold on;
    h1=plot(t,Q,'k','linewidth',1.5);
    for i = 1:size(flow_section,1)
        h2=plot(t(flow_section(i,1):flow_section(i,2)),...
            Q(flow_section(i,1):flow_section(i,2)),'r-','linewidth',1.5);
    end
    h3=plot(t,minimum_peak.*ones(size(t)),'k--');
    title('Selected rising limbs')
    legend([h1 h2 h3],{'Full flow series', 'Selected rising limbs', 'Minimum peak threshold'})
    % datetick('x')
    xlabel('Date')
    ylabel('Flow [mm/timestep]')
    fig_handles.RisingLimbs = fig;
end

end