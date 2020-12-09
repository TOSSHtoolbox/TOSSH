function [results] = calc_All(Q_mat, t_mat, P_mat, PET_mat, T_mat)
%calc_All calculates all signatures in the toolbox.
%   If a signature function can calculate multiple signatures
%   (e.g. sig_x_percentile) only one signature is calculated (e.g. Q95).
%   Note: This function is primarily intended to test all signatures.
%
%   INPUT
%   Q_mat: streamflow [mm/timestep] matrix (cell array)
%   t_mat: time [Matlab datenum] matrix (cell array)
%   P_mat: precipitation [mm/timestep] matrix (cell array)
%   PET_mat: pot. evapotranspiration [mm/timestep] matrix (cell array)
%   T_mat: temperature [degC] matrix (cell array)
%
%   OUTPUT
%   results: struc array with all results (each signature for each time
%       series and associated error strings)
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   % create consistent cell arrays
%   Q_mat = {data.Q};
%   t_mat = {data.t};
%   P_mat = {data.P};
%   PET_mat = {data.PET};
%   T_mat = {data.T};
%   results = calc_All(Q_mat,t_mat,P_mat,PET_mat,T_mat);
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 5
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true; % to be able to use t for time and T for temperature

% required input arguments
% Please input time series as a cell array of the following format:
% {x_1; x_2; ...; x_n}, where each entry (1, 2, ..., n) corresponds to one
% time series, e.g. from one catchment. For one catchment only, please
% input {x}. Example: {Q_1; Q_2; ...; Q_n} for streamflow.
addRequired(ip, 'Q_mat', @(Q_mat) iscell(Q_mat))
addRequired(ip, 't_mat', @(t_mat) iscell(t_mat))
addRequired(ip, 'P_mat', @(P_mat) iscell(P_mat))
addRequired(ip, 'PET_mat', @(PET_mat) iscell(PET_mat))
addRequired(ip, 'T_mat', @(T_mat) iscell(T_mat))

parse(ip, Q_mat, t_mat, P_mat, PET_mat, T_mat)

% calculate signatures

% initialise arrays
AC1 = NaN(size(Q_mat,1),1);
AC1_error_str = strings(size(Q_mat,1),1);
BaseflowRecessionK = NaN(size(Q_mat,1),1);
BaseflowRecessionK_error_str = strings(size(Q_mat,1),1);
BaseflowMagnitude = NaN(size(Q_mat,1),1);
BaseflowMagnitude_error_str = strings(size(Q_mat,1),1);
BFI = NaN(size(Q_mat,1),1);
BFI_error_str = strings(size(Q_mat,1),1);
EventGraphThresholds = NaN(size(Q_mat,1),10);
EventGraphThresholds_error_str = strings(size(Q_mat,1),1);
EventRR = NaN(size(Q_mat,1),1);
EventRR_error_str = strings(size(Q_mat,1),1);
FDC = cell(size(Q_mat,1),1);
FDC_error_str = strings(size(Q_mat,1),1);
FDC_slope = NaN(size(Q_mat,1),1);
FDC_slope_error_str = strings(size(Q_mat,1),1);
FlashinessIndex = NaN(size(Q_mat,1),1);
FlashinessIndex_error_str = strings(size(Q_mat,1),1);
HFD_mean = NaN(size(Q_mat,1),1);
HFD_mean_error_str = strings(size(Q_mat,1),1);
HFI_mean = NaN(size(Q_mat,1),1);
HFI_mean_error_str = strings(size(Q_mat,1),1);
MRC_SlopeChanges = cell(size(Q_mat,1),2);
MRC_SlopeChanges_error_str = strings(size(Q_mat,1),1);
PeakDistribution = NaN(size(Q_mat,1),1);
PeakDistribution_error_str = strings(size(Q_mat,1),1);
PQ_Curve = NaN(size(Q_mat,1),4);
PQ_Curve_error_str = strings(size(Q_mat,1),1);
Q_CoV = NaN(size(Q_mat,1),1);
Q_CoV_error_str = strings(size(Q_mat,1),1);
Q_mean = NaN(size(Q_mat,1),1);
Q_mean_error_str = strings(size(Q_mat,1),1);
Q_mean_monthly = NaN(size(Q_mat,1),1);
Q_mean_monthly_error_str = strings(size(Q_mat,1),1);
Q_7_day_max = NaN(size(Q_mat,1),1);
Q_7_day_max_error_str = strings(size(Q_mat,1),1);
Q_7_day_min = NaN(size(Q_mat,1),1);
Q_7_day_min_error_str = strings(size(Q_mat,1),1);
Q_skew = NaN(size(Q_mat,1),1);
Q_skew_error_str = strings(size(Q_mat,1),1);
Q_var = NaN(size(Q_mat,1),1);
Q_var_error_str = strings(size(Q_mat,1),1);
QP_elasticity = NaN(size(Q_mat,1),1);
QP_elasticity_error_str = strings(size(Q_mat,1),1);
RecessionParameters = NaN(size(Q_mat,1),2);
RecessionParameters_error_str = strings(size(Q_mat,1),1);
RecessionK_early = NaN(size(Q_mat,1),1);
RecessionK_early_error_str = strings(size(Q_mat,1),1);
Spearmans_rho = NaN(size(Q_mat,1),1);
Spearmans_rho_error_str = strings(size(Q_mat,1),1);
ResponseTime = NaN(size(Q_mat,1),1);
ResponseTime_error_str = strings(size(Q_mat,1),1);
RLD = NaN(size(Q_mat,1),1);
RLD_error_str = strings(size(Q_mat,1),1);
RR_Seasonality = NaN(size(Q_mat,1),1);
RR_Seasonality_error_str = strings(size(Q_mat,1),1);
SeasonalTranslation = NaN(size(Q_mat,1),1);
SeasonalTranslation_error_str = strings(size(Q_mat,1),1);
Recession_a_Seasonality = NaN(size(Q_mat,1),1);
Recession_a_Seasonality_error_str = strings(size(Q_mat,1),1);
SnowDayRatio = NaN(size(Q_mat,1),1);
SnowDayRatio_error_str = strings(size(Q_mat,1),1);
SnowStorage = NaN(size(Q_mat,1),1);
SnowStorage_error_str = strings(size(Q_mat,1),1);
StorageFraction = NaN(size(Q_mat,1),3);
StorageFraction_error_str = strings(size(Q_mat,1),1);
StorageFromBaseflow = NaN(size(Q_mat,1),1);
StorageFromBaseflow_error_str = strings(size(Q_mat,1),1);
TotalRR = NaN(size(Q_mat,1),1);
TotalRR_error_str = strings(size(Q_mat,1),1);
VariabilityIndex = NaN(size(Q_mat,1),1);
VariabilityIndex_error_str = strings(size(Q_mat,1),1);
Q95 = NaN(size(Q_mat,1),1);
Q95_error_str = strings(size(Q_mat,1),1);
high_Q_duration = NaN(size(Q_mat,1),1);
high_Q_duration_error_str = strings(size(Q_mat,1),1);
high_Q_frequency = NaN(size(Q_mat,1),1);
high_Q_frequency_error_str = strings(size(Q_mat,1),1);

% warning('off','all')

% loop over all catchments
for i = 1:size(Q_mat,1)
    
    if mod(i,1) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,size(Q_mat,1))
    end
    
    [AC1(i),~,AC1_error_str(i)] = sig_Autocorrelation(Q_mat{i},t_mat{i});
    [BaseflowRecessionK(i),~,BaseflowRecessionK_error_str(i)] = ...
        sig_BaseflowRecessionK(Q_mat{i},t_mat{i},'eps',median(Q_mat{i},'omitnan'));
    [BaseflowMagnitude(i),~,BaseflowMagnitude_error_str(i)] = sig_BaseflowMagnitude(Q_mat{i},t_mat{i});
    [BFI(i),~,BFI_error_str(i)] = sig_BFI(Q_mat{i},t_mat{i});
    [EventGraphThresholds(i,1),EventGraphThresholds(i,2),...
        EventGraphThresholds(i,3),EventGraphThresholds(i,4),...
        EventGraphThresholds(i,5),EventGraphThresholds(i,6),...
        EventGraphThresholds(i,7),EventGraphThresholds(i,8),...
        EventGraphThresholds(i,9),EventGraphThresholds(i,10),...
        ~,EventGraphThresholds_error_str(i)] = ...
        sig_EventGraphThresholds(Q_mat{i},t_mat{i},P_mat{i});
    [EventRR(i),~,EventRR_error_str(i)] = sig_EventRR(Q_mat{i},t_mat{i},P_mat{i});
    [FDC{i}(:,1), FDC{i}(:,2),~,FDC_error_str(i)] = sig_FDC(Q_mat{i},t_mat{i});
    [FDC_slope(i),~,FDC_slope_error_str(i)] = sig_FDC_slope(Q_mat{i},t_mat{i});
    [FlashinessIndex(i),~,FlashinessIndex_error_str(i)] = sig_FlashinessIndex(Q_mat{i},t_mat{i});
    [HFD_mean(i),~,HFD_mean_error_str(i)] = sig_HFD_mean(Q_mat{i},t_mat{i});
    [HFI_mean(i),~,HFI_mean_error_str(i)] = sig_HFI_mean(Q_mat{i},t_mat{i});
    [MRC_SlopeChanges{i,1},MRC_SlopeChanges{i,2},~,MRC_SlopeChanges_error_str(i)] = ...
        sig_MRC_SlopeChanges(Q_mat{i},t_mat{i});
    [PeakDistribution(i),~,PeakDistribution_error_str(i)] = sig_PeakDistribution(Q_mat{i},t_mat{i});
    [PQ_Curve(i,1),PQ_Curve(i,2),PQ_Curve(i,3),PQ_Curve(i,4),~,PQ_Curve_error_str(i)] = ...
        sig_PQ_Curve(Q_mat{i},t_mat{i},P_mat{i});
    [Q_CoV(i),~,Q_CoV_error_str(i)] = sig_Q_CoV(Q_mat{i},t_mat{i});
    [Q_mean(i),~,Q_mean_error_str(i)] = sig_Q_mean(Q_mat{i},t_mat{i});
    [Q_mean_monthly(i),~,Q_mean_monthly_error_str(i)] = sig_Q_mean_monthly(Q_mat{i},t_mat{i},1);
    [Q_7_day_max(i),~,Q_7_day_max_error_str(i)] = sig_Q_n_day_max(Q_mat{i},t_mat{i},7);
    [Q_7_day_min(i),~,Q_7_day_min_error_str(i)] = sig_Q_n_day_min(Q_mat{i},t_mat{i},7);
    [Q_skew(i),~,Q_skew_error_str(i)] = sig_Q_skew(Q_mat{i},t_mat{i});
    [Q_var(i),~,Q_var_error_str(i)] = sig_Q_var(Q_mat{i},t_mat{i});
    [QP_elasticity(i),~,QP_elasticity_error_str(i)] = sig_QP_elasticity(Q_mat{i},t_mat{i},P_mat{i});
    [RecessionParameters(i,:),~,~,RecessionParameters_error_str(i)] = ...
        sig_RecessionAnalysis(Q_mat{i},t_mat{i},'fit_individual',false);
    [RecessionK_early(i),~,RecessionK_early_error_str(i)] = sig_RecessionParts(Q_mat{i},t_mat{i},'early');
    [Spearmans_rho(i),~,Spearmans_rho_error_str(i)] = sig_RecessionUniqueness(Q_mat{i},t_mat{i});
    [ResponseTime(i),~,ResponseTime_error_str(i)] = sig_ResponseTime(Q_mat{i},t_mat{i},P_mat{i});
    [RLD(i),~,RLD_error_str(i)] = sig_RisingLimbDensity(Q_mat{i},t_mat{i});
    [RR_Seasonality(i),~,RR_Seasonality_error_str(i)] = sig_RR_Seasonality(Q_mat{i},t_mat{i},P_mat{i});
    [SeasonalTranslation(i,1),SeasonalTranslation(i,2),~,SeasonalTranslation_error_str(i)] = ...
        sig_SeasonalTranslation(Q_mat{i},t_mat{i},P_mat{i},PET_mat{i});
    [Recession_a_Seasonality(i),~,Recession_a_Seasonality_error_str(i)] = sig_SeasonalVarRecessions(Q_mat{i},t_mat{i});
    [SnowDayRatio(i),~,SnowDayRatio_error_str(i)] = sig_SnowDayRatio(Q_mat{i},t_mat{i},P_mat{i},T_mat{i});
    [SnowStorage(i),~,SnowStorage_error_str(i)] = sig_SnowStorage(Q_mat{i},t_mat{i},P_mat{i});
    [StorageFraction(i,1),StorageFraction(i,2),StorageFraction(i,3),~,StorageFraction_error_str(i)] = ... 
        sig_StorageFraction(Q_mat{i},t_mat{i},P_mat{i},PET_mat{i});
    [StorageFromBaseflow(i),~,StorageFromBaseflow_error_str(i)] = ...
        sig_StorageFromBaseflow(Q_mat{i},t_mat{i},P_mat{i},PET_mat{i});
    [TotalRR(i),~,TotalRR_error_str(i)] = sig_TotalRR(Q_mat{i},t_mat{i},P_mat{i});
    [VariabilityIndex(i),~,VariabilityIndex_error_str(i)] = sig_VariabilityIndex(Q_mat{i},t_mat{i});
    [Q95(i),~,Q95_error_str(i)] = sig_x_percentile(Q_mat{i},t_mat{i},95);
    [high_Q_duration(i),~,high_Q_duration_error_str(i)] = sig_x_Q_duration(Q_mat{i},t_mat{i},'high');
    [high_Q_frequency(i),~,high_Q_frequency_error_str(i)] = sig_x_Q_frequency(Q_mat{i},t_mat{i},'high');
    
end

% warning('on','all')

% add results to struct array
results.AC1 = AC1;
results.AC1_error_str = AC1_error_str;
results.BaseflowRecessionK = BaseflowRecessionK;
results.BaseflowRecessionK_error_str = BaseflowRecessionK_error_str;
results.BaseflowMagnitude = BaseflowMagnitude;
results.BaseflowMagnitude_error_str = BaseflowMagnitude_error_str;
results.BFI = BFI;
results.BFI_error_str = BFI_error_str;
results.EventGraphThresholds = EventGraphThresholds;
results.EventGraphThresholds_error_str = EventGraphThresholds_error_str;
results.EventRR = EventRR;
results.EventRR_error_str = EventRR_error_str;
results.FDC = FDC;
results.FDC_error_str = FDC_error_str;
results.FDC_slope = FDC_slope;
results.FDC_slope_error_str = FDC_slope_error_str;
results.FlashinessIndex = FlashinessIndex;
results.FlashinessIndex_error_str = FlashinessIndex_error_str;
results.HFD_mean = HFD_mean;
results.HFD_mean_error_str = HFD_mean_error_str;
results.HFI_mean = HFI_mean;
results.HFI_mean_error_str = HFI_mean_error_str;
results.MRC_SlopeChanges = MRC_SlopeChanges;
results.MRC_SlopeChanges_error_str = MRC_SlopeChanges_error_str;
results.PeakDistribution = PeakDistribution;
results.PeakDistribution_error_str = PeakDistribution_error_str;
results.PQ_Curve = PQ_Curve;
results.PQ_Curve_error_str = PQ_Curve_error_str;
results.Q_CoV = Q_CoV;
results.Q_CoV_error_str = Q_CoV_error_str;
results.Q_mean = Q_mean;
results.Q_mean_error_str = Q_mean_error_str;
results.Q_mean_monthly = Q_mean_monthly;
results.Q_mean_monthly_error_str = Q_mean_monthly_error_str;
results.Q_7_day_max = Q_7_day_max;
results.Q_7_day_max_error_str = Q_7_day_max_error_str;
results.Q_7_day_min = Q_7_day_min;
results.Q_7_day_min_error_str = Q_7_day_min_error_str;
results.Q_skew = Q_skew;
results.Q_skew_error_str = Q_skew_error_str;
results.Q_var = Q_var;
results.Q_var_error_str = Q_var_error_str;
results.QP_elasticity = QP_elasticity;
results.QP_elasticity_error_str = QP_elasticity_error_str;
results.RecessionParameters = RecessionParameters;
results.RecessionParameters_error_str = RecessionParameters_error_str;
results.RecessionK_early = RecessionK_early;
results.RecessionK_early_error_str = RecessionK_early_error_str;
results.Spearmans_rho = Spearmans_rho;
results.Spearmans_rho_error_str = Spearmans_rho_error_str;
results.ResponseTime = ResponseTime;
results.ResponseTime_error_str = ResponseTime_error_str;
results.RLD = RLD;
results.RLD_error_str = RLD_error_str;
results.RR_Seasonality = RR_Seasonality;
results.RR_Seasonality_error_str = RR_Seasonality_error_str;
results.SeasonalTranslation = SeasonalTranslation;
results.SeasonalTranslation_error_str = SeasonalTranslation_error_str;
results.Recession_a_Seasonality = Recession_a_Seasonality;
results.Recession_a_Seasonality_error_str = Recession_a_Seasonality_error_str;
results.SnowDayRatio = SnowDayRatio;
results.SnowDayRatio_error_str = SnowDayRatio_error_str;
results.SnowStorage = SnowStorage;
results.SnowStorage_error_str = SnowStorage_error_str;
results.StorageFraction = StorageFraction;
results.StorageFraction_error_str = StorageFraction_error_str;
results.StorageFromBaseflow = StorageFromBaseflow;
results.StorageFromBaseflow_error_str = StorageFromBaseflow_error_str;
results.TotalRR = TotalRR;
results.TotalRR_error_str = TotalRR_error_str;
results.VariabilityIndex = VariabilityIndex;
results.VariabilityIndex_error_str = VariabilityIndex_error_str;
results.Q95 = Q95;
results.Q95_error_str = Q95_error_str;
results.high_Q_duration = high_Q_duration;
results.high_Q_duration_error_str = high_Q_duration_error_str;
results.high_Q_frequency = high_Q_frequency;
results.high_Q_frequency_error_str = high_Q_frequency_error_str;

end
