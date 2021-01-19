function [P, PET, Q, T] = loadCatchmentCAMELSGB(ID,path)
%loadCatchmentCAMELSGB Loads hydro-meteorological time series (P, PET, Q, 
%   T) for CAMELS-GB format.
%
%   INPUT
%   ID: catchment ID
%   path: file path
%
%   OUTPUT
%   P: precipitation [mm/d]
%   PET: potential evapotranspiration (without interception correction) [mm/d]
%   Q: streamflow [mm/d]
%   T: T [degC]
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

file_ID = strcat(path,'CAMELS_GB_hydromet_timeseries_',num2str(ID),'_19701001-20150930.csv');

% date	precipitation	pet	temperature	discharge_spec	discharge_vol	peti	humidity	shortwave_rad	longwave_rad	windspeed
data = readtable(file_ID);

date = datenum(data.date);
Q_temp = data.discharge_spec;
P_temp = data.precipitation;
PET_temp = data.pet;
% PET_temp = data.peti; % with interception correction
T_temp = data.temperature;

Q = [date Q_temp];
P = [date P_temp];
PET = [date PET_temp];
T = [date T_temp];

end

