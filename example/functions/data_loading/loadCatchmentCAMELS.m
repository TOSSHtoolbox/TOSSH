function [P, PET, Q, T] = loadCatchmentCAMELS(ID,path_mod,path_obs,area)
%loadCatchmentCAMELS Loads hydro-meteorological time series (P, PET, Q, T).
%   We use the modelled time series (typically from October 1980 to
%   December 2014) to load P, PET, and T data (P and T are copied into
%   these files, but are observed values). We use the observed flow series
%   to load Q, which need to be converted to mm/day and adjusted as they
%   differ in length compared to the other time series.
%
%   INPUT
%   ID: catchment ID
%   path_mod: file path used to load P, PET, and T time series
%   path_obs: file path used to load flow series
%   area: catchment area (Gages)
%
%   OUTPUT
%   P: precipitation [mm/d]
%   PET: potential evapotranspiration (adjusted using standard coefficient
%   of 1.26) [mm/d]
%   Q: streamflow [mm/d]
%   T: T [°C]
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 4
    error('Not enough input arguments.')
end

foundID = false; % search through all folders until catchment is found
i_str_list = ["01","02","03","04","05","06","07","08","09","10","11",...
    "12","13","14","15","16","17","18"];
i = 0;
while foundID == false
    i = i+1;
    i_str = i_str_list(i);
    try
        % modelled time series
        file_ID_model = strcat(path_mod,i_str,'\',num2str(ID,'%08d'),'_05_model_output.txt');
        txt_data=fileread(file_ID_model);
        % model parameters
        file_ID_parameters = strcat(path_mod,i_str,'\',num2str(ID,'%08d'),'_05_model_parameters.txt');
        txt_para=fileread(file_ID_parameters);
        % observed time series
        file_ID_flow = strcat(path_obs,i_str,'\',num2str(ID,'%08d'),'_streamflow_qc.txt');
        txt_flow=fileread(file_ID_flow);
        foundID = true;
    catch
        disp('')
        foundID = false;
    end
end

% NOTE: NOW version 05 of Newman dataset
% YR MNTH DY HR SWE PRCP RAIM TAIR PET ET MOD_RUN OBS_RUN
data_model_cell = textscan(txt_data,...
    '%f %f %f %f %f %f %f %f %f %f %f %f', ...
    'Delimiter', '\t', 'HeaderLines', 1);

data_parameter_cell = textscan(txt_para,...
    '%s %f', ...
    'Delimiter', '\t', 'HeaderLines', 0);
PET_coefficient = data_parameter_cell{2}(41);

% GAGEID Year Month Day Streamflow(cubic feet per second) QC_flag
data_flow_cell = textscan(txt_flow,...
    '%f %f %f %f %f %s', ...
    'Delimiter', '\t', 'HeaderLines', 0);

Y = data_model_cell{1};
M = data_model_cell{2};
D = data_model_cell{3};
date = datenum(Y,M,D);

% Q_temp = data_model_cell{12};
Q_temp = data_flow_cell{5};
% cubicft/s to mm/day: Q = (q/35.3146667)*(86.4/area)
Q_temp(Q_temp==-999) = NaN;
Q_temp = (Q_temp./35.3146667).*(86.4/area);
% remove values before start of P, PET, and T time series 
Y_Q = data_flow_cell{2};
M_Q = data_flow_cell{3};
D_Q = data_flow_cell{4};
date_Q = datenum(Y_Q,M_Q,D_Q);
Q_temp(date_Q<date(1)) = []; 
% add NaNs at the end if Q time series is too short
Q_temp2 = Q_temp;
Q_temp = NaN(size(date));
Q_temp(1:length(Q_temp2)) = Q_temp2;

P_temp = data_model_cell{6};
PET_temp = data_model_cell{9}; 
T_temp = data_model_cell{8};

Q = [date Q_temp];
P = [date P_temp];
% PET = [date PET_temp]; % calibrated PET
PET = [date (1.26/PET_coefficient).*PET_temp]; % adjusted PET
T = [date T_temp];

end

