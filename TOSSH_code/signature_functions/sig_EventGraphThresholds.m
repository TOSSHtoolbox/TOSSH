function [IE_effect, SE_effect, IE_thresh_signif, IE_thresh, ...
    SE_thresh_signif, SE_thresh, SE_slope, ...
    Storage_thresh, Storage_thresh_signif, min_Qf_perc, ...
    error_flag, error_str, fig_handles] = sig_EventGraphThresholds(Q,t,P,varargin)
%%   sig_EventGraphThresholds calculates a variety of signatures related to saturation and infiltration excess flow.
%   Calculates a variety of signatures related to saturation excess (SE)
%   and infiltration excess (IE) flow, based on the relationship between
%   precipitation or antecedent conditions and flow response (volume or
%   peak magnitude) during storm events.
%
%   Notes:
%   Uses util_EventSeparation to separate individual storm events.
%   Enables plotting of the relationship between rainfall/flow quantities,
%   and calculation of regression coefficients and threshold strength.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   OPTIONAL
%   (All these parameters enable the user to tweak the definition of storm
%   events:)
%   min_termination: minimum termination time between storm events [hours],
%       default = 8
%   min_duration: minimum duration of storm [hours], default = 5
%   min_intensity_hour: minimum hourly rainfall [mm/hour], default = 2
%   min_intensity_day: minimum daily rainfall [mm/day], default = 10
%   min_intensity_hour_during: minimum rainfall allowed during an event
%       without contributing to termination time [mm/hour], default = 0.2
%   min_intensity_day_during: minimum rainfall allowed during an event
%       without contributing to termination time [mm/day], default = 1
%   max_recessiondays: maximum length of recession period after rainfall
%       [days] deemed to be part of the preceeding event (shorter if
%       another event starts), default = 1
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   IE_effect: infiltration excess (IE) importance [-]
%   SE_effect: saturation excess (SE) importance [-]
%   IE and SE importance based on their coefficients in regression
%   equations to predict event flow characteristics, adapted from
%   qualitative description in Estrany et al. (2010).
%
%   IE_thresh_signif: infiltration excess threshold significance [-]
%   IE_thresh: infiltration excess threshold location [mm/timestep]
%   Significance (using likelihood ratio test) and location of a threshold
%   in a plot of quickflow volume vs. maximum intensity, signifying IE
%   process, adapted from qualitative description in Ali et al. (2013);
%   IE is indicated when ie_thresh_sig < 0.05.
%
%   SE_thresh_signif: saturation excess threshold significance [-]
%   SE_thresh: saturation excess threshold location [mm]
%   SE_slope: saturation excess threshold above-threshold slope [mm/mm]
%   Significance, location, and above-threshold slope of a threshold in a
%   plot of quickflow volume vs. total precipitation. No threshold
%   indicates flow generation from riparian areas (Tani, 1997). Slope above
%   threshold indicates rate at which saturated areas expand (Tani, 1997;
%   Becker and McDonnell, 1998); SE is indicated when se_thresh_sig < 0.05.
%
%   Storage_thresh: storage/saturation excess threshold significance [-]
%   Storage_thresh_signif: storage/saturation excess threshold location [mm]
%   Significance and location of a threshold in a plot of quickflow volume
%   vs. antecedent precipitation index + total precipitation, adapted from
%   qualitative description in Ali et al. (2013) and McGrath et al. (2007);
%   SE is indicated when storage_thresh_sig < 0.05.
%
%   min_Qf_perc: minimum quickflow as a percentage of precipitation [%],
%   Indicates impermeable area contribution (qualitative description in
%   Becker and McDonnell, 1998).
%
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
%   [IE_effect, SE_effect] = sig_EventGraphThresholds(Q,t,P);
%   [IE_effect, SE_effect, IE_thresh_signif, IE_thresh, ...
%   SE_thresh_signif, SE_thresh, SE_slope, Storage_thresh, ...
%   Storage_thresh_signif, min_Qf_perc] = ...
%   sig_EventGraphThresholds(Q,t,P,'plot_results',true);
%
%   References
%   Ali, G., Oswald, C.J., Spence, C., Cammeraat, E.L., McGuire, K.J.,
%   Meixner, T. and Reaney, S.M., 2013. Towards a unified threshold-based
%   hydrological theory: necessary components and recurring challenges.
%   Hydrological Processes, 27(2), pp.313-318.
%   Becker, A. and McDonnell, J.J., 1998. Topographical and ecological
%   controls of runoff generation and lateral flows in mountain catchments.
%   IAHS Publications-Series of Proceedings and Reports-Intern Assoc
%   Hydrological Sciences, 248, pp.199-206.
%   Estrany, J., Garcia, C. and Batalla, R.J., 2010. Hydrological response
%   of a small mediterranean agricultural catchment. Journal of Hydrology,
%   380(1-2), pp.180-190.
%   McGrath, G.S., Hinz, C. and Sivapalan, M., 2007. Temporal dynamics of
%   hydrological threshold events. Hydrology and Earth System Sciences,
%   11(2), pp.923-938.
%   Mosley, M.P., 1979. Streamflow generation in a forested watershed, New
%   Zealand. Water Resources Research, 15(4), pp.795-806.
%   Tani, M., 1997. Runoff generation processes estimated from hydrological
%   observations on a steep forested hillslope with a thin soil layer.
%   Journal of Hydrology, 200(1-4), pp.84-109.
%   Wrede, S., Fenicia, F., Martinez-Carreras, N., Juilleret, J., Hissler,
%   C., Krein, A., Savenije, H.H., Uhlenbrook, S., Kavetski, D. and
%   Pfister, L., 2015. Towards more systematic perceptual model
%   development: a case study using 3 Luxembourgish catchments.
%   Hydrological Processes, 29(12), pp.2731-2750.
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
% flow time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) ((isnumeric(t) || isdatetime(t))) && ((size(t,1)==1 || size(t,2)==1)))
% precipitation time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'P', @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1))

% optional input arguments
addParameter(ip, 'min_termination', 24, @isnumeric) % minimum termination time
% (time between storms) in hours
addParameter(ip, 'min_duration', 5, @isnumeric) % minimum duration of storm in hour
addParameter(ip, 'min_intensity_hour', 2, @isnumeric) % minimum hourly rainfall (per hour)
addParameter(ip, 'min_intensity_day', 10, @isnumeric) % minimum daily rainfall (per day)
addParameter(ip, 'min_intensity_hour_during', 0.2, @isnumeric) % minimum timestep
% intensity allowed during storm event without contributing to termination time
addParameter(ip, 'min_intensity_day_during', 1, @isnumeric) % minimum timestep
% intensity allowed during storm event without contributing to termination time
addParameter(ip, 'max_recessiondays', 1, @isnumeric) % maximum number of days to allow recession after rain
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (1 graph)

parse(ip, Q, t, P, varargin{:})

min_termination = ip.Results.min_termination;
min_duration = ip.Results.min_duration;
min_intensity_hour = ip.Results.min_intensity_hour;
min_intensity_day = ip.Results.min_intensity_day;
min_intensity_hour_during = ip.Results.min_intensity_hour_during;
min_intensity_day_during = ip.Results.min_intensity_day_during;
max_recessiondays = ip.Results.max_recessiondays;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P);
if error_flag == 2
    IE_effect = NaN;
    SE_effect = NaN;
    IE_thresh_signif = NaN;
    IE_thresh = NaN;
    SE_thresh_signif = NaN;
    SE_thresh = NaN;
    SE_slope = NaN;
    Storage_thresh = NaN;
    Storage_thresh_signif = NaN;
    min_Qf_perc = NaN;
    return
end
timestep_factor = 1/days(timestep); % adjust for timestep
timestep = hours(timestep);

%% separate P and Q series into events

% run event separation algorithm
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[stormarray, error_flag, error_str, fig_handles] = util_EventSeparation(...
    datenum(t), P, timestep, min_termination, min_duration, ...
    min_intensity_hour, min_intensity_day, ...
    min_intensity_hour_during, min_intensity_day_during, ...
    max_recessiondays, plot_results);
if error_flag == 3
    IE_effect = NaN;
    SE_effect = NaN;
    IE_thresh_signif = NaN;
    IE_thresh = NaN;
    SE_thresh_signif = NaN;
    SE_thresh = NaN;
    SE_slope = NaN;
    Storage_thresh = NaN;
    Storage_thresh_signif = NaN;
    min_Qf_perc = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

% remove events where there is missing flow data
storm_array_check = true(size(stormarray,1),1);
for i = 1:size(stormarray,1)
    storm_array_check(i) = ~isnan(sum(Q(stormarray(i,1):stormarray(i,3))));
end
stormarray = stormarray(storm_array_check,:);

%% separate flow series into baseflow and quickflow
% create baseflow array
filter_parameter = exp(log(0.925)/timestep_factor);
B = util_LyneHollickFilter(Q, 'filter_parameter', filter_parameter, 'nr_passes', 3);

%% calculate P and Q characteristics for each storm
% cycle through storms and calculate the required storm characteristics
event_array = zeros(size(stormarray,1),9);
for i = 1:size(stormarray,1)
    % total event precipitation [mm] (Estrany names this as X3 coefficient)
    event_array(i,1) = sum(P(stormarray(i,1):stormarray(i,3)),'omitnan');
    
    % event average precipitation intensity [mm/timestep] (Estrany X4)
    event_array(i,2) = mean(P(stormarray(i,1):stormarray(i,3)),'omitnan');
    
    % event maximum precipitation intensity [mm/timestep] (Estrany X5; but they used 5 min data)
    event_array(i,3) = max(P(stormarray(i,1):stormarray(i,3)));
    
    % antecedent precipitation 3 days before [mm] (Estrany X1)
    time_antecedent = max(stormarray(i,1) - 3*24/timestep,1);
    event_array(i,4) = sum(P(time_antecedent:stormarray(i,1)),'omitnan');
    
    % antecedent precipitation 7 days before [mm] (Estrany X2)
    time_antecedent = max(stormarray(i,1) - 7*24/timestep,1);
    event_array(i,5) = sum(P(time_antecedent:stormarray(i,1)),'omitnan');
    
    % total event flow volume [mm]
    event_array(i,6) = sum(Q(stormarray(i,1):stormarray(i,3)),'omitnan');
    
    % total event quickflow volume [mm]
    event_array(i,7) = sum(Q(stormarray(i,1):stormarray(i,3)),'omitnan')-...
        sum(B(stormarray(i,1):stormarray(i,3)),'omitnan');
    
    % maximum event runoff [mm/timestep]
    event_array(i,8) = max(Q(stormarray(i,1):stormarray(i,3)));
    
    % antecedent precipitation index (API) using Mosley (1979) method
    % (sum Pi/i for 30 days)
    % get total antecedent time
    time_antecedent = max(stormarray(i,1) - 30*24/timestep,1);
    % form into daily blocks
    numdays = floor((stormarray(i,1) - time_antecedent)/(24/timestep));
    time_array = reshape(P(stormarray(i,1)-[1:1:(numdays*24/timestep)]),24/timestep,[]);
    time_array = sum(time_array,1);
    event_array(i,9) = sum(time_array./[1:length(time_array)],'omitnan').';
end

% calculate Becker and McDonnell (1998) signatures, minimum value of
% quickflow as a percentage of precipitation
min_Qf_perc = 100*(min(event_array(i,7)./event_array(i,1)));

% check if time series only consists of baseflow
if all(event_array(:,6)==0) || all(event_array(:,7)==0) || all(event_array(:,8)==0)
    IE_effect = NaN;
    SE_effect = NaN;
    IE_thresh_signif = NaN;
    IE_thresh = NaN;
    SE_thresh_signif = NaN;
    SE_thresh = NaN;
    SE_slope = NaN;
    Storage_thresh = NaN;
    Storage_thresh_signif = NaN;
    min_Qf_perc = NaN;
    error_flag = 3;
    error_str = ['Error: Event flow or quickflow volumes are zero for all events. ', error_str];
    return
end

%% calculate importance of storage vs. intensity characteristics
% this is done via multiple regressions to predict max flow and quickflow
% volumes (following Estrany et al., 2010)

% regression for maximum flow
flow_events = event_array(:,8) > 0;
if numel(flow_events(flow_events==1)) <= 1
    IE_effect = NaN;
    SE_effect = NaN;
    IE_thresh_signif = NaN;
    IE_thresh = NaN;
    SE_thresh_signif = NaN;
    SE_thresh = NaN;
    SE_slope = NaN;
    Storage_thresh = NaN;
    Storage_thresh_signif = NaN;
    min_Qf_perc = NaN;
    error_flag = 3;
    error_str = ['Error: Fewer than two suitable storm events found, try relaxing storm event criteria. ', error_str];
    return
end
mdl_max = stepwiselm(zscore(event_array(flow_events,1:5)),zscore(log(event_array(flow_events,8))),'Upper','linear','Verbose',0);
max_coeffs = mdl_max.Coefficients.Estimate;
max_vars = mdl_max.VariableInfo.InModel;

% regression for quickflow volume
flow_events = event_array(:,7) > 0;
if numel(flow_events(flow_events==1)) <= 1
    IE_effect = NaN;
    SE_effect = NaN;
    IE_thresh_signif = NaN;
    IE_thresh = NaN;
    SE_thresh_signif = NaN;
    SE_thresh = NaN;
    SE_slope = NaN;
    Storage_thresh = NaN;
    Storage_thresh_signif = NaN;
    min_Qf_perc = NaN;
    error_flag = 3;
    error_str = ['Error: Fewer than two suitable storm events found, try relaxing storm event criteria. ', error_str];
    return
end
mdl_qf = stepwiselm(zscore(event_array(flow_events,1:5)),zscore(log(event_array(flow_events,7))),'Upper','linear','Verbose',0);
qf_coeffs = mdl_qf.Coefficients.Estimate;
qf_vars = mdl_qf.VariableInfo.InModel;

% average scores for infiltration explanatory variables
ie_indices_max = and(max_vars(1:5),[0;1;1;0;0]);
ie_var_numbers_max = max_vars(1:5).*cumsum(max_vars(1:5));
ie_coeff_numbers_max = ie_var_numbers_max(ie_indices_max)+1;

ie_indices_qf = and(qf_vars(1:5),[0;1;1;0;0]);
ie_var_numbers_qf = qf_vars(1:5).*cumsum(qf_vars(1:5));
ie_coeff_numbers_qf = ie_var_numbers_qf(ie_indices_qf)+1;

IE_effect = mean([sum(max_coeffs(ie_coeff_numbers_max));sum(qf_coeffs(ie_coeff_numbers_qf))]);

% average scores for storage-related explanatory variables
se_indices_max = and(max_vars(1:5),[1;0;0;1;1]);
se_var_numbers_max = max_vars(1:5).*cumsum(max_vars(1:5));
se_coeff_numbers_max = se_var_numbers_max(se_indices_max)+1;

se_indices_qf = and(qf_vars(1:5),[1;0;0;1;1]);
se_var_numbers_qf = qf_vars(1:5).*cumsum(qf_vars(1:5));
se_coeff_numbers_qf = se_var_numbers_qf(se_indices_qf)+1;

SE_effect = mean([sum(max_coeffs(se_coeff_numbers_max));sum(qf_coeffs(se_coeff_numbers_qf))]);

%% check for thresholds and return locations and significance scores
% total precip against quickflow
[thresh_tp_qf,slope_tp_qf,slope_linear_tp_qf,p_value_tp_qf] = ...
    util_Threshold(event_array(:,1), event_array(:,7));
SE_thresh_signif = p_value_tp_qf;
SE_thresh = thresh_tp_qf;
SE_slope = slope_tp_qf;

% maximum intensity against quickflow
[thresh_mi_qf,slope_mi_qf,slope_linear_mi_qf,p_value_mi_qf] = ...
    util_Threshold(event_array(:,3), event_array(:,7));
IE_thresh_signif = p_value_mi_qf;
IE_thresh = thresh_mi_qf;

% API + total precipitation (proxy for storage) against quickflow
[thresh_st_qf,slope_st_qf,slope_linear_st_qf,p_value_st_qf] = ...
    util_Threshold(event_array(:,1)+event_array(:,9), event_array(:,7));
Storage_thresh_signif = p_value_st_qf;
Storage_thresh = thresh_st_qf;

%% optional plotting
if plot_results
    
    % plot flow, baseflow and events
    fig_events = figure('Position',[100 100 700 250]);
    title(['Flow, Baseflow and Event Separation'])
    dates_dt = t;%datetime(dates,'ConvertFrom','datenum');
    P_max = max(Q);
    hold on
    for i = 1:size(stormarray,1)
        p3=fill([dates_dt(stormarray(i,1)),dates_dt(stormarray(i,1)),...
            dates_dt(stormarray(i,2)),dates_dt(stormarray(i,2))],...
            [0, P_max, P_max, 0],'b','LineStyle','none');
        p4=fill([dates_dt(stormarray(i,2)),dates_dt(stormarray(i,2)),...
            dates_dt(stormarray(i,3)),dates_dt(stormarray(i,3))],...
            [0, P_max, P_max, 0],'c','LineStyle','none');
    end
    p1 = plot(t,Q,'r-');
    p2 = plot(t,B,'g-');
    p5=plot(dates_dt(:),P(:)./50,'k-','linewidth',1);

    xlabel('Time')
    ylabel('Flow [mm]')
    legend([p1, p2, p3, p4, p5],{'Flow','Baseflow','Events','Recessions','Rainfall/50'},'location','best')

    % plot total precip against quickflow, with threshold
    fig = figure('Position',[100 100 900 250]);
    subplot(1,3,1)
    title(['Quickflow vs P, p = ',num2str(round(p_value_tp_qf,6))])
    hold on
    p1 = plot(event_array(:,1),event_array(:,7),'bo');
    p2 = plot([min(event_array(:,1)),max(event_array(:,1))],...
        [min(event_array(:,1)),max(event_array(:,1))].*slope_linear_tp_qf,'k-');
    plot([0 thresh_tp_qf],[0 0],'g-')
    p3 = plot([thresh_tp_qf max(event_array(:,1))],...
        [0 (max(event_array(:,1))-thresh_tp_qf).*slope_tp_qf],'g-');
    xlabel('Event precipitation [mm]')
    ylabel('Event quickflow volume [mm]')
    legend([p1, p2, p3],{'Data','No threshold','With threshold'},'location','best')
    
    % maximum intensity against quickflow, with threshold
    subplot(1,3,2)
    title(['Quickflow vs intensity, p = ',num2str(round(p_value_mi_qf,6))])
    hold on
    p1 = plot(event_array(:,3),event_array(:,7),'bo');
    p2 = plot([min(event_array(:,3)),max(event_array(:,3))],...
        [min(event_array(:,3)),max(event_array(:,3))].*slope_linear_mi_qf,'k-');
    plot([0 thresh_mi_qf],[0 0],'g-')
    p3 = plot([thresh_mi_qf max(event_array(:,3))],...
        [0 (max(event_array(:,3))-thresh_mi_qf).*slope_mi_qf],'g-');
    xlabel('Event maximum intensity [mm/timestep]')
    ylabel('Event quickflow volume [mm]')
    legend([p1, p2, p3],{'Data','No threshold','With threshold'},'location','best')
    
    % API + total precipitation against quickflow, with threshold
    subplot(1,3,3)
    title(['Quickflow vs API + P, p = ',num2str(round(p_value_st_qf,6))])
    hold on
    x_val = event_array(:,1)+event_array(:,9);
    p1 = plot(x_val,event_array(:,7),'bo');
    p2 = plot([min(x_val),max(x_val)],[min(x_val),max(x_val)].*slope_linear_st_qf,'k-');
    plot([0 thresh_st_qf],[0 0],'g-')
    p3 = plot([thresh_st_qf max(x_val)],[0 (max(x_val)-thresh_st_qf).*slope_st_qf],'g-');
    xlabel('Event API + precipitation [mm]')
    ylabel('Event quickflow volume [mm]')
    legend([p1, p2, p3],{'Data','No threshold','With threshold'},'location','best')
    fig_handles.EventGraphThresholds = fig;
    
    % Wrede et al. (2015) say "high intensity storms that don't produce
    % flow imply no IE processes". We plot max intensity against quickflow
    % on a seasonal basis for the user to evaluate this.
    fig2 = figure('Position',[100 100 350 300]);
    % split events by season
    tmonth = month(t(stormarray(:,1)));
    i_djf = ismember(tmonth,[12,1,2]);
    i_mam = ismember(tmonth,[3,4,5]);
    i_jja = ismember(tmonth,[6,7,8]);
    i_son = ismember(tmonth,[9,10,11]);
    % plot per season
    plot(event_array(i_djf,3),event_array(i_djf,7),'bo')
    hold on
    plot(event_array(i_mam,3),event_array(i_mam,7),'gx')
    plot(event_array(i_jja,3),event_array(i_jja,7),'r^')
    plot(event_array(i_son,3),event_array(i_son,7),'ks')
    xlabel('Event maximum intensity [mm/timestep]')
    ylabel('Event quickflow volume [mm]')
    legend('DJF','MAM','JJA','SON','location','best')
    title('Quickflow vs intensity per season') % Wrede 2015
    fig_handles.EventGraphThresholdsSeasons = fig2;
    
end
