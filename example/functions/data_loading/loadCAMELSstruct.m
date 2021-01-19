function [CAMELS_data] = loadCAMELSstruct()
%loadCAMELSstruct Creates struct file with CAMELS data.
%   - Loads hydro-meteorological time series and catchment attributes
%   - Adjusts PET data to use standard Priestley-Taylor coefficient of 1.26
%   - Uses local paths (which contain large CAMELS files)
%   - Data can be found at: https://ral.ucar.edu/solutions/products/camels
%
%   INPUT
%
%   OUTPUT
%   CAMELS_data: struct file with CAMELS data
%
%   References
%   Newman, A.J., Clark, M.P., Sampson, K., Wood, A., Hay, L.E., Bock, A., 
%   Viger, R.J., Blodgett, D., Brekke, L., Arnold, J.R. and Hopson, T., 
%   2015. Development of a large-sample watershed-scale hydrometeorological 
%   data set for the contiguous USA: data set characteristics and 
%   assessment of regional variability in hydrologic model performance. 
%   Hydrology and Earth System Sciences, 19(1), p.209.
%   Addor, N., Newman, A.J., Mizukami, N. and Clark, M.P., 2017. The CAMELS
%   data set: catchment attributes and meteorology for large-sample
%   studies. Hydrology and Earth System Sciences (HESS), 21(10), 
%   pp.5293-5313.
%   N. Addor, A. Newman, M. Mizukami, and M. P. Clark, 2017. Catchment 
%   attributes for large-sample studies. Boulder, CO: UCAR/NCAR. 
%   https://doi.org/10.5065/D6G73C3Q
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

%% Specify paths
% We have to be in the TOSSH directory and the CAMELS data should be stored
% in a folder named CAMELS_v2.0 in the example_data directory. The  
% following folders are required:
% ./example/example_data/CAMELS_v2.0/camels_attributes_v2.0/camels_*.txt 
% (7 files; contain catchment attributes)
% ./example/example_data/CAMELS_v2.0/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/*/*_model_output.txt 
% (18 folders with >1000 files; contain forcing time series)
% ./example/example_data/CAMELS_v2.0/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/*/*_streamflow_qc.txt  
% (18 folders with 671 files; contain streamflow time series)

path_catchment_attributes = "./example/example_data/CAMELS_v2.0/camels_attributes_v2.0/camels_attributes_v2.0/";
path_modelled_time_series = "./example/example_data/CAMELS_v2.0/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/"; 
path_observed_time_series = "./example/example_data/CAMELS_v2.0/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/"; 

if ~(exist(path_catchment_attributes) == 7)
    error('Cannot find local path. You can download CAMELS from https://ral.ucar.edu/solutions/products/camels.')
elseif ~(exist(path_modelled_time_series) == 7)
    error('Cannot find local path. You can download CAMELS from https://ral.ucar.edu/solutions/products/camels.')
elseif ~(exist(path_observed_time_series) == 7)
    error('Cannot find local path. You can download CAMELS from https://ral.ucar.edu/solutions/products/camels.')
end

%% Load catchment attributes
% We first load the catchment attribute data which are saved in several txt 
% files.

% topography
% gauge_id;gauge_lat;gauge_lon;elev_mean;slope_mean;area_gages2;area_geospa_fabric
% file_ID_topo = fopen(strcat(path_catchment_attributes,'camels_topo.txt'),'r');
% file_topo = fread(file_ID_topo,'*char');
% file_topo = strrep(file_topo','NA','NaN'); % replace NA with NaN
% camels_topo_data = textscan(file_topo,'%f %f %f %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_topo);
camels_topo_data = readtable(strcat(path_catchment_attributes,'camels_topo.txt'));

% climate
% gauge_id;p_mean;pet_mean;p_seasonality;frac_snow;aridity;high_prec_freq;high_prec_dur;high_prec_timing;low_prec_freq;low_prec_dur;low_prec_timing
% file_ID_clim = fopen(strcat(path_catchment_attributes,'camels_clim.txt'),'r');
% file_clim = fread(file_ID_clim,'*char');
% file_clim = strrep(file_clim','NA','NaN'); % replace NA with NaN
% camels_climate_data = textscan(file_clim,'%f %f %f %f %f %f %f %f %s %f %f %s',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_clim);
camels_climate_data = readtable(strcat(path_catchment_attributes,'camels_clim.txt'));

% hydrology
% gauge_id;q_mean;runoff_ratio;slope_fdc;baseflow_index;stream_elas;q5;q95;high_q_freq;high_q_dur;low_q_freq;low_q_dur;zero_q_freq;hfd_mean
% file_ID_hydro = fopen(strcat(path_catchment_attributes,'camels_hydro.txt'),'r');
% file_hydro = fread(file_ID_hydro,'*char');
% file_hydro = strrep(file_hydro','NA','NaN'); % replace NA with NaN
% camels_hydro_data = textscan(file_hydro,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_hydro);
camels_hydro_data = readtable(strcat(path_catchment_attributes,'camels_hydro.txt'));

% soils
% gauge_id;soil_depth_pelletier;soil_depth_statsgo;soil_porosity;soil_conductivity;max_water_content;sand_frac;silt_frac;clay_frac;water_frac;organic_frac;other_frac
% file_ID_soil = fopen(strcat(path_catchment_attributes,'camels_soil.txt'),'r');
% file_soil = fread(file_ID_soil,'*char');
% file_soil = strrep(file_soil','NA','NaN'); % replace NA with NaN
% camels_soil_data = textscan(file_soil,'%f %f %f %f %f %f %f %f %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_soil);
camels_soil_data = readtable(strcat(path_catchment_attributes,'camels_soil.txt'));

% geology
% gauge_id;geol_1st_class;glim_1st_class_frac;geol_2nd_class;glim_2nd_class_frac;carbonate_rocks_frac;geol_porostiy;geol_permeability
% file_ID_geol = fopen(strcat(path_catchment_attributes,'camels_geol.txt'),'r');
% file_geol = fread(file_ID_geol,'*char');
% file_geol = strrep(file_geol','NA','NaN'); % replace NA with NaN
% camels_geol_data = textscan(file_geol,'%f %s %f %s %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_geol);
camels_geol_data = readtable(strcat(path_catchment_attributes,'camels_geol.txt'));

% vegetation
% gauge_id;frac_forest;lai_max;lai_diff;gvf_max;gvf_diff;dom_land_cover_frac;dom_land_cover;root_depth_50;root_depth_99
% file_ID_vege = fopen(strcat(path_catchment_attributes,'camels_vege.txt'),'r');
% file_vege = fread(file_ID_vege,'*char');
% file_vege = strrep(file_vege','NA','NaN'); % replace NA with NaN
% camels_vege_data = textscan(file_vege,'%f %f %f %f %f %f %f %s %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_vege);
camels_vege_data = readtable(strcat(path_catchment_attributes,'camels_vege.txt'));

% We then create arrays containing the catchment attributes and metadata.
gauge_id = camels_climate_data{:,1};

% topography
gauge_lat = camels_topo_data.gauge_lat;
gauge_lon = camels_topo_data.gauge_lon;
elev_mean = camels_topo_data.elev_mean;
slope_mean = camels_topo_data.slope_mean;
area_gages2 = camels_topo_data.area_gages2;
area_geospa_fabric = camels_topo_data.area_geospa_fabric;

% climate
p_mean = camels_climate_data.p_mean;
pet_mean = camels_climate_data.pet_mean;
p_seasonality = camels_climate_data.p_seasonality;
frac_snow = camels_climate_data.frac_snow;
aridity = camels_climate_data.aridity;
high_prec_freq = camels_climate_data.high_prec_freq;
high_prec_dur = camels_climate_data.high_prec_dur;
high_prec_timing = camels_climate_data.high_prec_timing;
low_prec_freq = camels_climate_data.low_prec_freq;
low_prec_dur = camels_climate_data.low_prec_dur;
low_prec_timing = camels_climate_data.low_prec_timing;

% hydrology
q_mean = camels_hydro_data.q_mean;
runoff_ratio = camels_hydro_data.runoff_ratio;
slope_fdc = camels_hydro_data.slope_fdc;
baseflow_index = camels_hydro_data.baseflow_index;
stream_elas = camels_hydro_data.stream_elas;
q5 = camels_hydro_data.q5;
q95 = camels_hydro_data.q95;
high_q_freq = camels_hydro_data.high_q_freq;
high_q_dur = camels_hydro_data.high_q_dur;
low_q_freq = camels_hydro_data.low_q_freq;
low_q_dur = camels_hydro_data.low_q_dur;
zero_q_freq = camels_hydro_data.zero_q_freq;
hfd_mean = camels_hydro_data.hfd_mean;

% soil
soil_depth_pelletier = camels_soil_data.soil_depth_pelletier;
soil_depth_statsgo = camels_soil_data.soil_depth_statsgo;
soil_porosity = camels_soil_data.soil_porosity;
soil_conductivity = camels_soil_data.soil_conductivity;
max_water_content = camels_soil_data.max_water_content;
sand_frac = camels_soil_data.sand_frac;
silt_frac = camels_soil_data.silt_frac;
clay_frac = camels_soil_data.clay_frac;
water_frac = camels_soil_data.water_frac;
organic_frac = camels_soil_data.organic_frac;
other_frac = camels_soil_data.other_frac;

% geology
geol_1st_class = camels_geol_data.geol_1st_class;
glim_1st_class_frac = camels_geol_data.glim_1st_class_frac;
geol_2nd_class = camels_geol_data.geol_2nd_class;
glim_2nd_class_frac = camels_geol_data.glim_2nd_class_frac;
carbonate_rocks_frac = camels_geol_data.carbonate_rocks_frac;
geol_porosity = camels_geol_data.geol_porostiy;
geol_permeability = camels_geol_data.geol_permeability;

% vegetation
frac_forest = camels_vege_data.frac_forest;
lai_max = camels_vege_data.lai_max;
lai_diff = camels_vege_data.lai_diff;
gvf_max = camels_vege_data.gvf_max;
gvf_diff = camels_vege_data.gvf_diff;
dom_land_cover_frac = camels_vege_data.dom_land_cover_frac;
dom_land_cover = camels_vege_data.dom_land_cover;
root_depth_50 = camels_vege_data.root_depth_50;
root_depth_99 = camels_vege_data.root_depth_99;

%% Load hydro-meteorological time series
% To extract the time series, we loop over all catchments. We also
% calculate the completeness of the flow records.
flow_perc_complete = NaN(length(gauge_id),1); % completeness of flow record
P = cell(length(gauge_id),1); % precipitation
PET = cell(length(gauge_id),1); % potential evapotranspiration
Q = cell(length(gauge_id),1); % streamflow
T = cell(length(gauge_id),1); % temperature

fprintf('Loading catchment data...\n')
for i = 1:length(gauge_id) % loop over all catchments
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,length(gauge_id))
    end
    
    ID = gauge_id(i);
    [P{i}, PET{i}, Q{i}, T{i}] = loadCatchmentCAMELS(...
        ID,path_modelled_time_series,path_observed_time_series,area_gages2(i));
    flow_perc_complete(i) = 100*(1-sum(isnan(Q{i}(:,2)))./length(Q{i}(:,2)));
    
end

%% Create struct file with catchment attributes and time series
% We create a Matlab structure file containing all the data, i.e. the
% hydro-meteorological time series and the catchment attributes and
% metadata.

% topography
CAMELS_data.gauge_id = gauge_id;
CAMELS_data.gauge_lat = gauge_lat;
CAMELS_data.gauge_lon = gauge_lon;
CAMELS_data.elev_mean = elev_mean;
CAMELS_data.slope_mean = slope_mean;
CAMELS_data.area_gages2 = area_gages2;
CAMELS_data.area_geospa_fabric = area_geospa_fabric;

% climate
CAMELS_data.p_mean = p_mean;
CAMELS_data.pet_mean = pet_mean;
CAMELS_data.p_seasonality = p_seasonality;
CAMELS_data.frac_snow = frac_snow;
CAMELS_data.aridity = aridity;
CAMELS_data.high_prec_freq = high_prec_freq;
CAMELS_data.high_prec_dur = high_prec_dur;
CAMELS_data.high_prec_timing = high_prec_timing;
CAMELS_data.low_prec_freq = low_prec_freq;
CAMELS_data.low_prec_dur = low_prec_dur;
CAMELS_data.low_prec_timing = low_prec_timing;

% hydrology
CAMELS_data.q_mean = q_mean;
CAMELS_data.runoff_ratio = runoff_ratio;
CAMELS_data.slope_fdc = slope_fdc;
CAMELS_data.baseflow_index = baseflow_index;
CAMELS_data.stream_elas = stream_elas;
CAMELS_data.q5 = q5;
CAMELS_data.q95 = q95;
CAMELS_data.high_q_freq = high_q_freq;
CAMELS_data.high_q_dur = high_q_dur;
CAMELS_data.low_q_freq = low_q_freq;
CAMELS_data.low_q_dur = low_q_dur;
CAMELS_data.zero_q_freq = zero_q_freq;
CAMELS_data.hfd_mean = hfd_mean;

% soil
CAMELS_data.soil_depth_pelletier = soil_depth_pelletier;
CAMELS_data.soil_depth_statsgo = soil_depth_statsgo;
CAMELS_data.soil_porosity = soil_porosity;
CAMELS_data.soil_conductivity = soil_conductivity;
CAMELS_data.max_water_content = max_water_content;
CAMELS_data.sand_frac = sand_frac;
CAMELS_data.silt_frac = silt_frac;
CAMELS_data.clay_frac = clay_frac;
CAMELS_data.water_frac = water_frac;
CAMELS_data.organic_frac = organic_frac;
CAMELS_data.other_frac = other_frac;

% geology
CAMELS_data.geol_1st_class = geol_1st_class;
CAMELS_data.glim_1st_class_frac = glim_1st_class_frac;
CAMELS_data.geol_2nd_class = geol_2nd_class;
CAMELS_data.glim_2nd_class_frac = glim_2nd_class_frac;
CAMELS_data.carbonate_rocks_frac = carbonate_rocks_frac;
CAMELS_data.geol_porosity = geol_porosity;
CAMELS_data.geol_permeability = geol_permeability;

% vegetation
CAMELS_data.frac_forest = frac_forest;
CAMELS_data.lai_max = lai_max;
CAMELS_data.lai_diff = lai_diff;
CAMELS_data.gvf_max = gvf_max;
CAMELS_data.gvf_diff = gvf_diff;
CAMELS_data.dom_land_cover_frac = dom_land_cover_frac;
CAMELS_data.dom_land_cover = dom_land_cover;
CAMELS_data.root_depth_50 = root_depth_50;
CAMELS_data.root_depth_99 = root_depth_99;

% hydro-meteorological time series
CAMELS_data.flow_perc_complete = flow_perc_complete;
CAMELS_data.P = P;
CAMELS_data.PET = PET;
CAMELS_data.Q = Q;
CAMELS_data.T = T;

% To avoid loading CAMELS from the original files every time, we could also 
% save the struct file and just load that.
% save('./example/example_data/CAMELS_data.mat','CAMELS_data')

end