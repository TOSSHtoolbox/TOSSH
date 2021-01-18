%% TOSSH workflow 4 - calculation of signatures for CAMELS US catchments
%
%   This script shows how to use TOSSH to calculate the Addor et al. (2018)
%   signatures using the CAMELS dataset. Note that this workflow can be
%   slow and requires sufficient RAM.
%   
%   The CAMELS dataset can be
%   downloaded from https://ral.ucar.edu/solutions/products/camels 
%   and needs to be placed in the right directory.  
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
%   Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and 
%   Clark, M.P., 2018. A ranking of hydrological signatures based on their 
%   predictability in space. Water Resources Research, 54(11), 
%   pp.8792-8812.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

close all
% clear all
clc

%% Download and extract CAMELS data
% First, we need to download and extract the CAMELS data from:
% https://ral.ucar.edu/solutions/products/camels
% From the link above you can download six different folders. We need three 
% of them to run this workflow:
% basin_timeseries_v1p2_modelOutput_daymet (CAMELS time series meteorology, observed flow, meta data (.zip))
% basin_timeseries_v1p2_metForcing_obsFlow (CAMELS time series Daymet forced model output (.zip))
% camels_attributes_v2.0 (CAMELS Attributes (.zip))
% The data should be stored in a folder named CAMELS_v2.0 located in the 
% example_data directory of TOSSH (./example/example_data/CAMELS_v2.0). 

%% Add directories to path
% We navigate to the TOSSH directory and add it to the Matlab path. This is 
% important to ensure that we can work with relative paths. If we already 
% are in this directory, we can use the pwd command:
mydir = pwd;
% Alternatively, we can specify my_dir manually:
% mydir = 'D:/Sebastian/Documents/MATLAB/TOSSH';
cd(mydir)
addpath(genpath(mydir));

% We also specify the path where to save our results and figures:
results_path = './example/results/';
fig_path = './example/results/images';

%% Load CAMELS data
% We now load the CAMELS data into a struct file for easy handling with 
% Matlab. We have to be in the TOSSH directory and the CAMELS data should 
% be stored in a folder named CAMELS_v2.0 in the example_data directory.  

% The following folders are required:
% ./example/example_data/CAMELS_v2.0/camels_attributes_v2.0/camels_*.txt 
% (7 files; contain catchment attributes)
% ./example/example_data/CAMELS_v2.0/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/*/*_model_output.txt 
% (18 folders with >1000 files; contain forcing time series)
% ./example/example_data/CAMELS_v2.0/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2\/usgs_streamflow/*/*_streamflow_qc.txt  
% (18 folders with 671 files; contain streamflow time series)

% Loading the data might take a few minutes.
CAMELS_data = loadCAMELSstruct(); 
% Note that you can also save the struct file to avoid loading the data
% anew every time you want to work with them.

%% Calculate signatures for CAMELS catchments using TOSSH
% To use the calculation function calc_Addor.m, we need to create cell
% arrays containing the time series. We use cell arrays since some time
% series might have different lengths. While the length of each row in the 
% cell array can vary, the cell arrays containing the t, Q, P, and PET data
% need to have exactly the same dimensions. We first initialise the cell
% arrays.
n_CAMELS = length(CAMELS_data.gauge_id);
t_mat = cell(n_CAMELS,1);
Q_mat = cell(n_CAMELS,1);
P_mat = cell(n_CAMELS,1);
PET_mat = cell(n_CAMELS,1);

% We then loop over all catchments and extract the time series for each
% catchment. We crop each time series to be from October 1989 to September 
% 2009 (see Addor et al., 2017). 
fprintf('Creating data matrix...\n')
for i = 1:n_CAMELS 
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,n_CAMELS)
    end
    
    t = datetime(CAMELS_data.Q{i}(:,1),'ConvertFrom','datenum');
    Q = CAMELS_data.Q{i}(:,2);    
    P = CAMELS_data.P{i}(:,2);
    PET = CAMELS_data.PET{i}(:,2);
    
    % get subperiod, here from 1 Oct 1989 to 30 Sep 2009
    indices = 1:length(t); 
    start_ind = indices(t==datetime(1989,10,1));
    % in case time series starts after 1 Oct 1989
    if isempty(start_ind); start_ind = 1; end 
    end_ind = indices(t==datetime(2009,9,30));    
    t = t(start_ind:end_ind);
    Q = Q(start_ind:end_ind);
    P = P(start_ind:end_ind);
    PET = PET(start_ind:end_ind);
    % calculate completeness during sub-period
    flow_perc_complete = 100*(1-sum(isnan(Q))./length(Q));
    CAMELS_data.flow_perc_complete(i) = flow_perc_complete;
    
    t_mat{i} = t;
    Q_mat{i} = Q;
    P_mat{i} = P;
    PET_mat{i} = PET;
    
end

% We can now use the calculation function calc_Addor.m to calculate all the
% Addor et al. (2018) signatures.
fprintf('Calculating signatures...\n')
CAMELS_signatures = calc_Addor(Q_mat, t_mat, P_mat);
% Besides the signature values, the function also returns a list with 
% warnings and error messages. Some warnings are specific to a certain 
% signature. For example: "Warning: FDC slope could not be calculated, 
% probably because flow is intermittent." tells us why sometimes NaN is
% returned when we calculate FDC_slope. 
fprintf(CAMELS_signatures.FDC_slope_error_str(294)+"\n")
% Most warnings come from our data check and indicate that there are some
% NaN values in the time series. 
fprintf(CAMELS_signatures.FDC_slope_error_str(339)+"\n")

% We can save the results as mat file which can be easily loaded into
% Matlab. Alternatively, we can save the results as txt file.
save(strcat(results_path,'CAMELS_signatures.mat'),...
    '-struct','CAMELS_signatures')
writetable(struct2table(CAMELS_signatures),...
    strcat(results_path,'CAMELS_signatures.txt'))

%% Compare TOSSH signatures to CAMELS signatures
% We can compare the signatures contained in CAMELS with the signatures 
% calculated here to see if we get the same results.
makeScatterPlot(CAMELS_signatures,CAMELS_data,-10) 
saveFig(gcf,strcat('TOSSH_scatter_plot_US'),fig_path,'-dpdf') 
% The results show that for most signatures, our code matches the CAMELS
% data within the limits of small differences in signature definition. 
% The differences in the slope of the FDC are due to an error in the CAMELS 
% data (personal communication with Nans Addor).

%% Further analysis of resulting signatures
% We can also use the results to see how some signatures relate to 
% catchment attributes. For example, we can check how correlated the runoff
% ratio is with the aridity index and the snow fraction.
figure('pos',[100 100 400 300])
scatter(CAMELS_data.aridity,CAMELS_signatures.TotalRR,25,CAMELS_data.frac_snow,'filled')
caxis([0 0.8]); 
xlim([0 4])
c = colorbar;
title(c,'Snow fraction [-]')
xlabel('Aridity [-]')
ylabel('Runoff ratio [-]')

%% Further information 
% Further information can be found in the online documentation: 
% https://TOSSHtoolbox.github.io/TOSSH/ and in the other example scripts.
