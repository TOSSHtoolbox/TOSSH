function [Recession_Parameters, recession_month, error_flag, error_str, fig_handles] = ...
    sig_RecessionAnalysis(Q, t, varargin)
%sig_RecessionAnalysis calculates recession parameters.
%   dQ/dt = -a*Q^b
%   Fits power law function to recession segments and returns recession
%   parameters (see Brutsaert and Nieber, 1977; Roques et al., 2017; and
%   Jachens et al., 2020).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   recession_length: min. length of recession segments [days], default = 5
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default = 0
%   start_of_recession: define start of recession when baseflow filter
%       rejoins the curve ("baseflow"), or after hydrograph peak ("peak")
%   filter_par: smoothing parameter of Lyne-Hollick filter to determine
%      start of recession (higher = later recession start), default = 0.925
%   fit_individual: fit each individual recession segment
%   fitting_type: fit non-linear or linear curve ('nonlinear','linear')
%       reservoir), etc.
%   dQdt_method: method for dQ/dt calculation, default = 'ETS'
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   Recession_Parameters: matrix with parameters alpha, beta (=1 for 
%       exponential fit) for each recession segment.
%       Recession_Parameters(:,1): a, scaling parameter
%       Recession_Parameters(:,2): b, parameter of non-linearity
%   recession_month: approx. month of recession
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
%   [para_mat,~] = sig_RecessionAnalysis(Q,t);
%   [para_mat,recession_month] = sig_RecessionAnalysis(Q,t,...
%   'plot_results',true,'fit_individual',false,'fitting_type','linear');
%
%   References
%	Brutsaert, W. and Nieber, J.L., 1977. Regionalized drought flow
%   hydrographs from a mature glaciated plateau. Water Resources Research,
%   13(3), pp.637-643.
%   Roques, C., Rupp, D.E. and Selker, J.S., 2017. Improved streamflow
%   recession parameter estimation with attention to calculation of? dQ/dt.
%   Advances in Water Resources, 108, pp.29-43.
%   Jachens, E.R., Rupp, D.E., Roques, C. and Selker, J.S., 2020. Recession
%   analysis revisited: Impacts of climate on parameter estimation.
%   Hydrology and Earth System Sciences, 24(3), pp.1159-1170.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

%   NOTE: Post-processing of Recession Parameters
%
%   Post-processing of Recession parameters can be done as:
%   RecessionParameters_a(i) = median((RecessionParametersTemp(:,1)),'omitnan');
%   RecessionParameters_b(i) = median(RecessionParametersTemp(:,2),'omitnan');
%   RecessionParametersT0Temp = 1./(RecessionParametersTemp(:,1).*median(Q_mat{i}(Q_mat{i}>0),'omitnan').^(RecessionParametersTemp(:,2)-1));
%   ReasonableT0 = and(RecessionParametersTemp(:,2)>0.5,RecessionParametersTemp(:,2)<5);
%   RecessionParameters_T0(i) = median(RecessionParametersT0Temp(ReasonableT0),'omitnan');

%   RecessionParameters_T0 is characteristic timescale of recessions at median flow; It can be obtained 
%   by fitting a line to the dQ/dt versus Q point cloud in log-log space for each individual recession, 
%   with  Q scaled by median Q; T0 is the median value of −1/intercept (McMillan et al., 2021). 

%   RecessionParameters_T0 can be derived from RecessionParameters_a and RecessionParameters_b as follows (McMillan et al, 2014): 
%   Let Qhat = Q/Qmedian, then:
%       dQhat/dt = - (1/T0) * Qhat^b
%   Separating constant terms,
%       1/Qmedian * dQ/dt = - (1/T0) * Q^b * (1/Qmedian)^b
%       dQ/dt = - (1/T0) * Q^b * (1/Qmedian) ^ (b-1)
%   As our estimate of a and b parameters are from:
%       dQ/dt = - a * Q^b
%   It leads to: 
%       a = (1/T0) * (1/Qmedian) ^ (b-1)
%       T0 = (1/a) * (1/Qmedian) ^ (b-1)
%       T0 = 1 / (a * Qmedian ^ (b-1)) # This corresponds to the code for calculating RecessionParametersT0Temp
%
%   McMillan, H., Gueguen, M., Grimon, E., Woods, R., Clark, M., & Rupp, D. E. (2014). 
%   Spatial variability of hydrological processes and model structure diagnostics in a 50 km2 catchment. 
%   Hydrological Processes, 28(18), 4896-4913. https://doi.org/10.1002/hyp.9988


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
% flow section (amount of timesteps) to be declared a recession
addParameter(ip, 'n_start', 1, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'peak', @ischar) % defines start of a recession
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of 
% Lyne-Hollick Filter to determine start of recession (higher = later recession start)
addParameter(ip, 'fit_individual', true, @islogical) % fit individual recessions or point cloud
addParameter(ip, 'fitting_type', 'linear', @ischar) % nonlinear or linear fit
addParameter(ip, 'dQdt_method', 'ETS', @ischar) % how to calculate dQ/dt
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
n_start = ip.Results.n_start;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
filter_par = ip.Results.filter_par;
fit_individual = ip.Results.fit_individual;
fitting_type = ip.Results.fitting_type;
dQdt_method = ip.Results.dQdt_method;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    Recession_Parameters = NaN(1,2);
    recession_month = NaN;
    return
end

% calculate signature

% get recession segments
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[flow_section, error_flag, error_str, fig_handles] = util_RecessionSegments(Q, t, ...
    'recession_length', recession_length, 'eps', eps, ...
    'filter_par', filter_par, 'plot_results', plot_results, ...
    'start_of_recession', start_of_recession, 'n_start', n_start);
if error_flag == 3
    Recession_Parameters = NaN(1,2);
    recession_month = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

% get flow rate gradient and corresponding flows
[dQdt, Qm, flow_section, R2] = ...
    util_dQdt(Q, t, flow_section, 'method', dQdt_method);

% get recession month
date_tmp = datevec(t(floor(mean(flow_section,2))));
recession_month = date_tmp(:,2);

% calculate recession parameters
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
if ~fit_individual
    rec = ~isnan(Qm);
    [Recession_Parameters(1), Recession_Parameters(2), error_flag, error_str] = ...
        util_FitPowerLaw(Qm(rec), dQdt(rec), fitting_type, R2(rec));
    recession_month = NaN; % no recession month since we only fit a single curve
else
    Recession_Parameters = NaN(size(flow_section,1),2);
    for i = 1:size(flow_section,1)
        rec = [flow_section(i,1):flow_section(i,2)]'; % get recession
        goodrec = ~isnan(Qm(rec));
        [Recession_Parameters(i,1), Recession_Parameters(i,2), error_flag, error_str] = ...
            util_FitPowerLaw(Qm(rec(goodrec)), dQdt(rec(goodrec)), fitting_type, R2(rec(goodrec)));
    % remove recessions with unrealistic (negative) parameter values
    if Recession_Parameters(i,1) <= 0 || Recession_Parameters(i,2) <= 0
       Recession_Parameters(i,:) = NaN;
    end
    end
end
if error_flag == 3
    Recession_Parameters = NaN(1,2);
    recession_month = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

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
        rec = [flow_section(i,1):flow_section(i,2)]'; % get recession
        Q_tmp = Qm(rec);
        dQdt_tmp = dQdt(rec);
        if fit_individual
            % date_vec = datevec(t(rec));
            % ind = floor(mean(date_vec(:,2))); % get approx. month
            ind = recession_month(i);
            plot(Q_tmp,-dQdt_tmp,'.','color',colour_mat_seasons(ind,:),'linewidth',2)
            plot(Q_tmp,Recession_Parameters(i,1).*Q_tmp.^Recession_Parameters(i,2),'color',colour_mat_seasons(ind,:))
        else
            plot(Q_tmp,-dQdt_tmp,'b.','linewidth',2)
        end
    end
    if fit_individual
        legend([p1 p2 p3 p4],{'MAM','JJA','SON','DJF'},'box','off','Location','best');
    end
    
    if ~fit_individual
        rec = ~isnan(Qm);
        plot(sort(Qm(rec)),Recession_Parameters(1).*sort(Qm(rec)).^Recession_Parameters(2),...
            '-','color','k','linewidth',1) %,'DisplayName','Fit'
        str = (sprintf('-dQ/dt = %.2f Q^{%.1f} \n',Recession_Parameters(1),Recession_Parameters(2)));
        title(str);
    else
        [~, ind] = min(abs(Recession_Parameters(:,2) - median(Recession_Parameters(:,2),'omitnan')));
        str = (sprintf('median: -dQ/dt = %.2f Q^{%.1f}',Recession_Parameters(ind,1),Recession_Parameters(ind,2)));       
        title(str);
    end
    xlabel('Q [mm/timestep]')
    ylabel('-dQ/dt [mm/timestep^2]') 
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    fig_handles.RecessionAnalysis = fig;
end

end

