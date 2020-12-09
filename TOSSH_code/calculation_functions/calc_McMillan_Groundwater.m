function [results] = calc_McMillan_Groundwater(Q_mat, t_mat, P_mat, PET_mat, varargin)
%calc_McMillan_Groundwater calculates various groundwater signatures.
%   Calculates 15 signatures from McMillan (2020), related to groundwater 
%   storage, groundwater dynamics and baseflow. These signatures come from 
%   previous experimental studies that link watershed or hillslope 
%   processes to streamflow reponse dynamics. Some signatures are 
%   implemented direct from the original papers, others are interpreted
%   from a qualitative description in the paper.
%
%   INPUT
%   Q_mat: streamflow [mm/timestep] matrix (cell array)
%   t_mat: time [Matlab datenum] matrix (cell array)
%   P_mat: precipitation [mm/timestep] matrix (cell array)
%   PET_mat: pot. evapotranspiration [mm/timestep] matrix (cell array)
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
%   results = calc_McMillan_Groundwater(Q_mat,t_mat,P_mat,PET_mat);
%
%   References
%   McMillan, H., 2020. Linking hydrologic signatures to hydrologic
%   processes: A review. Hydrological Processes, 34(6), pp.1393-1409.
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
% Please input time series as a cell array of the following format:
% {x_1; x_2; ...; x_n}, where each entry (1, 2, ..., n) corresponds to one
% time series, e.g. from one catchment. For one catchment only, please
% input {x}. Example: {Q_1; Q_2; ...; Q_n} for streamflow.
addRequired(ip, 'Q_mat', @(Q_mat) iscell(Q_mat))
addRequired(ip, 't_mat', @(t_mat) iscell(t_mat))
addRequired(ip, 'P_mat', @(P_mat) iscell(P_mat))
addRequired(ip, 'PET_mat', @(PET_mat) iscell(PET_mat))

% optional input arguments
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q_mat, t_mat, P_mat, PET_mat, varargin{:})

plot_results = ip.Results.plot_results;

% calculate signatures

% initialise arrays

% Section: Groundwater 
% Signatures referring to double peaks in streamflow are not coded due to
% subjectivity in identification of the peaks.
% Total runoff ratio describes overall water loss to deep groundwater.
TotalRR = NaN(size(Q_mat,1),1);
TotalRR_error_str = strings(size(Q_mat,1),1);
% Ratio between summer and winter runoff ratios, Low ratios show high
% bedrock permeability, Pfister et al. (2017).
RR_Seasonality = NaN(size(Q_mat,1),1);
RR_Seasonality_error_str = strings(size(Q_mat,1),1);
% Low event runoff ratios show show rapid vertical drainage of water to
% groundwater, McMillan et al (2011), Noguchi et al. (1997).
EventRR = NaN(size(Q_mat,1),1);
EventRR_error_str = strings(size(Q_mat,1),1);
% Ratio of active to total storage volumes, low ratios show permeable
% bedrock and high total storage, Pfister et al. (2017).
StorageFraction = NaN(size(Q_mat,1),3);
StorageFraction_error_str = strings(size(Q_mat,1),1);

% Section: Storage (especially groundwater)
% Seasonal variations in recession rate (Shaw and Riha; 2012).
Recession_a_Seasonality = NaN(size(Q_mat,1),1);
Recession_a_Seasonality_error_str = strings(size(Q_mat,1),1);
% Average storage from average baseflow and storage-discharge relationship 
% (McNamara et al., 2011).
AverageStorage = NaN(size(Q_mat,1),1);
AverageStorage_error_str = strings(size(Q_mat,1),1);
% Recession analysis parameters approximate storage-discharge relationship.
% This fits a single relationship to the point cloud.
RecessionParameters = NaN(size(Q_mat,1),2);
RecessionParameters_error_str = strings(size(Q_mat,1),1);
% Changes of slope in a master recession curve or recession analysis plot
% indicate multiple linear reservoirs.
MRC_num_segments = NaN(size(Q_mat,1),1);
MRC_num_segments_error_str = strings(size(Q_mat,1),1);
% First, steep section of the master recession curve = storage that is
% quickly depleted (Estrany et al., 2010).
First_Recession_Slope = NaN(size(Q_mat,1),1);
% Mid section of the master recession curve = Water retention capacity of 
% the watershed (Estrany et al., 2010).
Mid_Recession_Slope = NaN(size(Q_mat,1),1);
% Non-uniqueness in the recession relationship (McMillan et al., 2011;
% Harman et al., 2009). Tested using spearman rank correlation on Q vs dQdt.
Spearmans_rho = NaN(size(Q_mat,1),1);
Spearmans_rho_error_str = strings(size(Q_mat,1),1);
% Ratio between event and total Runoff coefficient: Low = Large storage 
% capacity (Blume et al., 2008).
EventRR_TotalRR_ratio = NaN(size(Q_mat,1),1);
% Variability index of flow. Low variability index shows higher water 
% storage (Estrany et al., 2010).
VariabilityIndex = NaN(size(Q_mat,1),1);
VariabilityIndex_error_str = strings(size(Q_mat,1),1);

% Section: Baseflow 
% Visual inspection of hydrograph for stable baseflow: not coded.
% Baseflow index indicates baseflow proportion and baseflow residence time 
% (Yilmaz et al., 2008; Bulygina et al., 2009; Hrachowitz et al., 2014).
BFI = NaN(size(Q_mat,1),1);
BFI_error_str = strings(size(Q_mat,1),1);
% Baseflow recession constant K, slower recessions show greater groundwater
% influence and longer subsurface flow paths (Safeeq et al., 2013).
BaseflowRecessionK = NaN(size(Q_mat,1),1);
BaseflowRecessionK_error_str = strings(size(Q_mat,1),1);

% loop over all catchments
for i = 1:size(Q_mat,1)
    
    % Section: Groundwater
    [TotalRR(i),~,TotalRR_error_str(i)] = sig_TotalRR(Q_mat{i},t_mat{i},P_mat{i});
    [RR_Seasonality(i),~,RR_Seasonality_error_str(i)] = sig_RR_Seasonality(Q_mat{i}, t_mat{i}, P_mat{i});
    [EventRR(i),~,EventRR_error_str(i)] = sig_EventRR(Q_mat{i},t_mat{i},P_mat{i});
    [StorageFraction(i,1),StorageFraction(i,2),StorageFraction(i,3),~,StorageFraction_error_str(i)] = ...
        sig_StorageFraction(Q_mat{i},t_mat{i},P_mat{i},PET_mat{i});
    
    % Section: Storage (especially groundwater)
    [Recession_a_Seasonality(i),~,Recession_a_Seasonality_error_str(i)] = ...
        sig_SeasonalVarRecessions(Q_mat{i},t_mat{i},'eps',median(Q_mat{i},'omitnan')*0.001,'plot_results',plot_results);
    [AverageStorage(i),~,AverageStorage_error_str(i)] = ...
        sig_StorageFromBaseflow(Q_mat{i},t_mat{i},P_mat{i},PET_mat{i},'plot_results',plot_results);
    [RecessionParameters(i,:),~,~,RecessionParameters_error_str(i)] = ...
        sig_RecessionAnalysis(Q_mat{i},t_mat{i},'fit_individual',false,'plot_results',plot_results);
    [MRC_num_segments(i),Segment_slopes,~,MRC_num_segments_error_str(i)] = ...
        sig_MRC_SlopeChanges(Q_mat{i},t_mat{i},'plot_results',plot_results,'eps',0.001*median(Q_mat{i},'omitnan'));
    First_Recession_Slope(i) = Segment_slopes(1);
    if length(Segment_slopes) >= 2
        Mid_Recession_Slope(i) = Segment_slopes(2);
    end
    [Spearmans_rho(i),~,Spearmans_rho_error_str(i)] = sig_RecessionUniqueness(Q_mat{i},t_mat{i});
    EventRR_TotalRR_ratio(i) = EventRR(i)/TotalRR(i);
    [VariabilityIndex(i),~,VariabilityIndex_error_str(i)] = sig_VariabilityIndex(Q_mat{i},t_mat{i});
    
    % Section: Baseflow
    [BFI(i),~,BFI_error_str(i)] = sig_BFI(Q_mat{i},t_mat{i},'method','UKIH','parameters',5);
    [BaseflowRecessionK(i),~,BaseflowRecessionK_error_str(i)] = ...
        sig_BaseflowRecessionK(Q_mat{i},t_mat{i},'eps',0.001*median(Q_mat{i},'omitnan')); 
    
end

% add results to struct array
results.TotalRR = TotalRR;
results.TotalRR_error_str = TotalRR_error_str;
results.EventRR = EventRR;
results.EventRR_error_str = EventRR_error_str;
results.RR_Seasonality = RR_Seasonality;
results.RR_Seasonality_error_str = RR_Seasonality_error_str;
results.StorageFraction = StorageFraction;
results.StorageFraction_error_str = StorageFraction_error_str;
results.Recession_a_Seasonality = Recession_a_Seasonality;
results.Recession_a_Seasonality_error_str = Recession_a_Seasonality_error_str;
results.AverageStorage = AverageStorage;
results.AverageStorage_error_str = AverageStorage_error_str;
results.RecessionParameters = RecessionParameters;
results.RecessionParameters_error_str = RecessionParameters_error_str;
results.MRC_num_segments = MRC_num_segments;
results.MRC_num_segments_error_str = MRC_num_segments_error_str;
results.BFI = BFI;
results.BFI_error_str = BFI_error_str;
results.BaseflowRecessionK = BaseflowRecessionK;
results.BaseflowRecessionK_error_str = BaseflowRecessionK_error_str;
results.First_Recession_Slope = First_Recession_Slope;
results.Mid_Recession_Slope = Mid_Recession_Slope;
results.Spearmans_rho = Spearmans_rho;
results.Spearmans_rho_error_str = Spearmans_rho_error_str;
results.EventRR_TotalRR_ratio = EventRR_TotalRR_ratio;
results.VariabilityIndex = VariabilityIndex;
results.VariabilityIndex_error_str = VariabilityIndex_error_str;

end
