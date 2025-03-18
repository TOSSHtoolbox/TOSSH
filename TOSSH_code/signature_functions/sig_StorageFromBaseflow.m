function [AverageStorage, error_flag, error_str, fig_handles] = ...
    sig_StorageFromBaseflow(Q, t, P, PET, varargin)
%sig_StorageFromBaseflow calculates average storage from average baseflow.
%   Uses a water balance approach to calculate daily changes in storage, 
%   then finds the relationship between storage and discharge, then
%   estimates average storage from average baseflow (see McNamara et al.,
%   2011 and Peters and Aulenbach, 2011).
%
%   Notes:
%   Modifies the method by conditioning AET on PET and soil moisture rather
%   than assuming AET = PET.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   PET: potential evapotranspiration [mm/timestep]
%   OPTIONAL
%   start_water_year: first month of water year, default = 10 (October)
%   field_capacity: field capacity [mm]
%   plot_results: whether to plot results, default = false
%   recession_length: min. length of recessions (days), default = 5
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default = 0
%       (note that large eps values can lead to problematic recession
%       selection)
%
%   OUTPUT
%   AverageStorage: storage [mm]
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
%   PET = data.PET;
%   AverageStorage = sig_StorageFromBaseflow(Q, t, P, PET);
%   AverageStorage = sig_StorageFromBaseflow(Q, t, P, PET, 'plot_results', true);
%
%   References
%   McNamara, J.P., Tetzlaff, D., Bishop, K., Soulsby, C., Seyfried, M., 
%   Peters, N.E., Aulenbach, B.T. and Hooper, R., 2011. Storage as a metric 
%   of catchment comparison. Hydrological Processes, 25(21), pp.3364-3371.
%   Peters, N.E. and Aulenbach, B.T., 2011. Water storage at the Panola 
%   Mountain research watershed, Georgia, USA. Hydrological Processes, 
%   25(25), pp.3878-3889.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 4
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
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'PET', @(PET) isnumeric(PET) && (size(PET,1)==1 || size(PET,2)==1))

% optional input arguments
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 1) && (x <= 12) && floor(x)==x;
addParameter(ip, 'start_water_year', 10, validationFcn) % when does the water year start? Default: 10
addParameter(ip, 'field_capacity', [], @isnumeric) % field capacity for scaling PET to AET
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)
addParameter(ip, 'recession_length', 5, @isnumeric) % length of decreasing
% flow in days to be declared a recession
addParameter(ip, 'n_start', 1, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period

parse(ip, Q, t, P, PET, varargin{:})
start_water_year = ip.Results.start_water_year;
field_capacity = ip.Results.field_capacity;
plot_results = ip.Results.plot_results;
recession_length = ip.Results.recession_length;
n_start = ip.Results.n_start;
eps = ip.Results.eps;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P, 'PET', PET);
if error_flag == 2
    AverageStorage = NaN;
    return
end

% calculate signature

% get rid of NaN values (temporarily)
isn = (isnan(Q) | isnan(P) | isnan(PET)); % store NaN indices
% replace NaN days with mean to have a roughly closed water balance
Q(isn) = mean(Q,'omitnan'); 
P(isn) = mean(P,'omitnan'); 
PET(isn) = mean(PET,'omitnan'); 

% estimate storage
[S, ~] = util_StorageAndAET(Q, t, P, PET, 'field_capacity', field_capacity);

% extract baseflow periods from the data
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[flow_section, error_flag, error_str] = util_RecessionSegments(Q, t, ...
    'recession_length', recession_length, 'n_start', n_start, 'eps', eps);
if error_flag == 3
    AverageStorage = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

num_points = sum(flow_section(:,2)-flow_section(:,1)+1);
storage_discharge = zeros(num_points,2);
storage_discharge_datetime = NaT(num_points,1);

counter = 0;

for i = 1:size(flow_section,1)
    storage_discharge(counter+1:counter+flow_section(i,2)-flow_section(i,1)+1,1) = ...
        S(flow_section(i,1):flow_section(i,2));
    storage_discharge(counter+1:counter+flow_section(i,2)-flow_section(i,1)+1,2) = ...
        Q(flow_section(i,1):flow_section(i,2));
    storage_discharge_datetime(counter+1:counter+flow_section(i,2)-flow_section(i,1)+1,1) = ...
        t(flow_section(i,1):flow_section(i,2));
    counter = counter + flow_section(i,2)-flow_section(i,1)+1;
end

good_points = find(storage_discharge(:,2)>0);
if numel(good_points)<size(storage_discharge,1)
    error_flag = 1;
    error_str = ['Warning: Ignoring zero discharge values in storage-discharge fit. ', error_str];
end
storage_discharge = storage_discharge(good_points,:);
storage_discharge_datetime = storage_discharge_datetime(good_points,:);

[water_year] = util_WaterYear(storage_discharge_datetime(:,1), 'WY_start_month', start_water_year);

% fit baseflow-storage relationship for combined water years, each with
% different intercept but a single slope
% S = m*ln(Q) + b
WY_unique = unique(water_year);
p_mat = zeros(length(storage_discharge),length(WY_unique));
for i = 1:length(storage_discharge)
    p_mat(i,find(WY_unique==water_year(i)))=1;
end
linFcn = @(p,x) p(1).*x + p_mat*p(2:end)';
p0 = [20, 10*ones(1,size(p_mat,2))];
fit_nonlin = fitnlm(log(storage_discharge(:,2)),storage_discharge(:,1),linFcn,p0);
% retrieve the slope and intercept fits
slope_fit = (fit_nonlin.Coefficients.Estimate(1));
intercept_coeffs = fit_nonlin.Coefficients.Estimate(2:end);

% adjust initial storage for each water year based on fit coefficients
revised_storage = storage_discharge(:,1);
for i = 1:length(WY_unique)
    % get all recession points for that water year and adjust initial storage
    wy_ind=find(water_year==WY_unique(i));
    revised_storage(wy_ind) = revised_storage(wy_ind)-intercept_coeffs(i);
end

% find average adjusted storage from average baseflow; average baseflow = 
% average of the 7-day minimum streamflow (McNamara et al., 2011)
baseflow = movmin(Q,7*(duration(24,0,0)/timestep));
av_baseflow = mean(baseflow);
AverageStorage = slope_fit.*log(av_baseflow) - min(revised_storage);

% optional plotting 
if plot_results
    % figure combining all the years 
    fig = figure('Position',[100 100 350 300]); hold on
    % scatter plot of baseflow vs adjusted storage
    gscatter(storage_discharge(:,2),revised_storage- min(revised_storage),water_year,[],[],10);
    xlabel('Baseflow [mm]')
    ylabel('Relative catchment storage [mm]')
    % fit and plot combined trendline
    fit_baseflow = sort(storage_discharge(:,2));
    fit_storage = slope_fit.*log(fit_baseflow);
    plot(fit_baseflow,fit_storage-min(revised_storage),'k-',...
        'DisplayName','Fitted line','Linewidth',1.5)
    legend('location','best')
    fig_handles.StorageFromBaseflow = fig;
end

end
