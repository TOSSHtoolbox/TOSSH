function [FDC, Q_sorted, error_flag, error_str, fig_handles] = sig_FDC(Q, t, varargin)
%sig_FDC calculates flow duration curve (FDC).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   FDC: Flow duration curve, i.e. exceedance probabilities [-]
%   Q_sorted: sorted streamflow values [mm/timestep]
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
%   [FDC, Q_sorted] = sig_FDC(Q,t,'plot_results',true);
%   
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1)) 
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1)) 

% optional input arguments
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, varargin{:})
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    FDC = NaN(size(Q));
    Q_sorted = NaN(size(Q));
    return
end

% calculate signature
% get ranks as a proxy for exceedance probabilities  
Q_tmp = Q(~isnan(Q)); % remove NaN values 
Q_sorted = sort(Q_tmp);
Q_ranked = [1:length(Q_tmp)]'; % give unique (random) rank to every measurement
FDC = 1 - Q_ranked./length(Q_ranked); % flow duration curve with unique ranks

% add warning for intermittent streams
if ~isempty(Q_tmp(Q_tmp==0))
    error_flag = 2;
    error_str = ['Warning: Flow is zero at least once (intermittent flow). ', error_str];
end

% optional plotting
if plot_results
    fig = figure('pos',[100 100 350 300]); hold on
    hold on    
    plot(FDC,Q_sorted,'linewidth',1.0)
    ylabel('Q [mm/timestep]')    
    xlabel('Exceedance probability [-]')
    set(gca,'YScale','log')
    fig_handles.FDC = fig;
end

end

