function [RLD, error_flag, error_str, fig_handles] = ...
    sig_RisingLimbDensity(Q, t, varargin)
%sig_RisingLimbDensity calculates rising limb density (RLD).
%   Calculates the rising limb density, the ratio between the number of
%   rising limbs and the total amount of timesteps the hydrograph is
%   rising (see e.g. Sawicz et al., 2011).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   rising_limb_length: length of increasing flow section (days) to be 
%       declared a rising limb, default = 1
%   eps: allowed increase in flow during rising limb, default = 0
%   minimum_peak: minimum peak to be counted as rising limb (peak size is
%       defined as difference between end and start of rising limb),
%       default = 0
%   plot_results: whether to plot results, default = 'false'
%
%   OUTPUT
%   RLD: rising limb density [1/timestep]
%   % rising_limb_month: approx. month of rising limb (not implemented)
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
%   RLD = sig_RisingLimbDensity(Q,t);
%   RLD = sig_RisingLimbDensity(Q,t,'plot_results',true,'rising_limb_length',2);
%
%   References
%   Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G.,
%   2011. Catchment classification: empirical analysis of hydrologic
%   similarity based on catchment function in the eastern USA. Hydrology
%   and Earth System Sciences, 15(9), pp.2895-2911.
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
addParameter(ip, 'rising_limb_length', 1, @isnumeric)% length of increasing flow section (timesteps) to be declared a rising limb
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

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    RLD = NaN;
    return
end

% calculate signature

% identify all rising limbs
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[flow_section, error_flag, error_str] = util_RisingLimbs(...
    Q, t, 'eps', eps, 'minimum_peak', minimum_peak, ...
    'rising_limb_length', rising_limb_length, 'plot_results', plot_results);
if error_flag == 3
    RLD = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

RLD = 1./mean(flow_section(:,2)-flow_section(:,1));

% % get rising limb month (not implemented)
% date_tmp = datevec(t(floor(mean(flow_section,2))));
% RLD_month = date_tmp(:,2);

% optional plotting
if plot_results
    fig = figure('Position',[100 100 350 300]); hold on
     colour_mat_seasons = [...
        0 0 1;  0 0 1;...
        0 1 0; 0 1 0; 0 1 0;...
        1 0 0; 1 0 0; 1 0 0;...
        1 1 0; 1 1 0; 1 1 0; ...
        0 0 1];
    p1=plot(0,0,'.','Color',[0 1 0]); 
    p2=plot(0,0,'.','Color',[1 0 0]);  
    p3=plot(0,0,'.','Color',[1 1 0]); 
    p4=plot(0,0,'.','Color',[0 0 1]);
    for i = 1:size(flow_section,1)
        RL = [flow_section(i,1):flow_section(i,2)]'; % get rising limb
        Q_tmp = Q(RL);
        t_tmp = t(RL);
        date_vec = datevec(t_tmp);
        ind = floor(median(date_vec(:,2))); % get approx. month
        plot(1:length(Q_tmp),Q_tmp,'-','color',colour_mat_seasons(ind,:),'linewidth',1)
    end
    xlabel('Time since start of rising limb [timestep]') 
    ylabel('Q [mm/timestep]') 
    set(gca,'YScale','log')
    legend([p1 p2 p3 p4],{'MAM','JJA','SON','DJF'},'box','off','Location','best');
    fig_handles.RisingLimbDensity = fig;
end

end
