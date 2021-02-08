function [results] = calc_Sawicz(Q_mat, t_mat, P_mat, T_mat, varargin)
%calc_Sawicz calculates signatures from Sawicz et al. (2011).
%   Sawicz et al. (2011) use 6 signatures drawn largely from Yadav et al. 
%   (2007), that are chosen to be uncorrelated and to be linked to 
%   catchment function. The signatures are used to analyze hydrological  
%   similarity between catchments, and link the resulting clusters to 
%   physical and climate attributes.
%
%   INPUT
%   Q_mat: streamflow [mm/timestep] matrix (cell array)
%   t_mat: time [Matlab datenum] matrix (cell array)
%   P_mat: precipitation [mm/timestep] matrix (cell array)
%   T_mat: temperature [degC] matrix (cell array)
%   OPTIONAL
%   start_water_year: first month of water year, default = 10 (October)
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
%   T_mat = {data.T}; 
%   results = calc_Sawicz(Q_mat,t_mat,P_mat,T_mat);
%
%   References
%	Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G.,
%   2011. Catchment classification: empirical analysis of hydrologic 
%   similarity based on catchment function in the eastern USA. Hydrology 
%   and Earth System Sciences, 15(9), pp.2895-2911.
%   Yadav, M., Wagener, T. and Gupta, H., 2007. Regionalization of 
%   constraints on expected watershed response behavior for improved
%   predictions in ungauged basins. Advances in Water Resources, 30(8), 
%   pp.1756-1774.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 4
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
addRequired(ip, 'T_mat', @(T_mat) iscell(T_mat))

% optional input arguments
addParameter(ip, 'start_water_year', 10, @isnumeric) % when does the water year start? Default: 10

parse(ip, Q_mat, t_mat, P_mat, T_mat, varargin{:})
start_water_year = ip.Results.start_water_year;

% initialise arrays
Total_RR = NaN(size(Q_mat,1),1);
Total_RR_error_str = strings(size(Q_mat,1),1);
FDC_slope = NaN(size(Q_mat,1),1);
FDC_slope_error_str = strings(size(Q_mat,1),1);
BFI = NaN(size(Q_mat,1),1);
BFI_error_str = strings(size(Q_mat,1),1);
QP_elasticity = NaN(size(Q_mat,1),1);
QP_elasticity_error_str = strings(size(Q_mat,1),1);
SnowDayRatio = NaN(size(Q_mat,1),1);
SnowDayRatio_error_str = strings(size(Q_mat,1),1);
RLD = NaN(size(Q_mat,1),1);
RLD_error_str = strings(size(Q_mat,1),1);

% loop over all catchments
for i = 1:size(Q_mat,1)  
    
    [Total_RR(i),~,Total_RR_error_str(i)] = sig_TotalRR(Q_mat{i},t_mat{i},P_mat{i});    
    [FDC_slope(i),~,FDC_slope_error_str(i)] = sig_FDC_slope(Q_mat{i},t_mat{i});    
    [BFI(i),~,BFI_error_str(i)] = sig_BFI(Q_mat{i},t_mat{i});    
    [QP_elasticity(i),~,QP_elasticity_error_str(i)] = sig_QP_elasticity(...
        Q_mat{i},t_mat{i},P_mat{i},'start_water_year',start_water_year);    
    [SnowDayRatio(i),~,SnowDayRatio_error_str(i)] = sig_SnowDayRatio(...
        Q_mat{i},t_mat{i},P_mat{i},T_mat{i});    
    [RLD(i),~,RLD_error_str(i)] = sig_RisingLimbDensity(Q_mat{i},t_mat{i});
    
end

% add results to struct array
results.Total_RR = Total_RR;
results.Total_RR_error_str = Total_RR_error_str;
results.FDC_slope = FDC_slope;
results.FDC_slope_error_str = FDC_slope_error_str;
results.BFI = BFI;
results.BFI_error_str = BFI_error_str;
results.QP_elasticity = QP_elasticity;
results.QP_elasticity_error_str = QP_elasticity_error_str;
results.SnowDayRatio = SnowDayRatio;
results.SnowDayRatio_error_str = SnowDayRatio_error_str;
results.RLD = RLD;
results.RLD_error_str = RLD_error_str;
    
end
