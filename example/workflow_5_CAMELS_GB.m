%% TOSSH workflow 5 - calculation of signatures for CAMELS GB catchments
%
%   This script shows how to use TOSSH to calculate various signatures
%   using the CAMELS GB dataset (Coxon et al., 2020). Note that this
%   workflow can be slow and requires sufficient RAM.
%
%   The CAMELS GB dataset can be downloaded from
%   https://catalogue.ceh.ac.uk/documents/8344e4f3-d2ea-44f5-8afa-86d2987543a9
%   and needs to be placed in the right directory.
%
%   References
%   Coxon, G., Addor, N., Bloomfield, J.P., Freer, J., Fry, M., Hannaford,
%   J., Howden, N.J., Lane, R., Lewis, M., Robinson, E.L. and Wagener, T.,
%   2020. CAMELS-GB: Hydrometeorological time series and landscape
%   attributes for 671 catchments in Great Britain. Earth System Science
%   Data Discussions, pp.1-34.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

close all
% clear all
clc

%% Download and extract CAMELS GB data
% First, we need to download and extract the CAMELS GB data from:
% https://catalogue.ceh.ac.uk/documents/8344e4f3-d2ea-44f5-8afa-86d2987543a9
% The data should be stored in a folder named CAMELS_GB located in the 
% example_data directory of TOSSH (.\example\example_data\CAMELS_GB). 

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

%% Load CAMELS GB data
% We now load the CAMELS GB data into a struct file for easy handling with
% Matlab. We have to be in the TOSSH directory and the CAMELS GB data  
% should be stored in a folder named CAMELS_GB in the example_data 
% directory.  

% The following folders are required:
% ./example/example_data/CAMELS_GB/data/CAMELS_GB_*_attributes.csv
% (8 files; contain catchment attributes)
% ./example/example_data/CAMELS_GB/data/timeseries/CAMELS_GB_hydromet*.csv 
% (671 files; contain forcing and streamflow time series) 

% Loading the datamight take a few minutes.
CAMELS_GB_data = loadCAMELSGBstruct();
% Note that you can also save the struct file to avoid loading the data
% anew every time you want to work with them.

%% Calculate signatures for CAMELS GB catchments using TOSSH
% To use the calculation function calc_Addor.m, we need to create cell
% arrays containing the time series. We use cell arrays since some time
% series might have different lengths. While the length of each row in the
% cell array can vary, the cell arrays containing the t, Q, P, and PET data
% need to have exactly the same dimensions. We first initialise the cell
% arrays.
n_CAMELS_GB = length(CAMELS_GB_data.gauge_id);
t_mat = cell(n_CAMELS_GB,1);
Q_mat = cell(n_CAMELS_GB,1);
P_mat = cell(n_CAMELS_GB,1);
PET_mat = cell(n_CAMELS_GB,1);

% We then loop over all catchments and extract the time series for each
% catchment.
fprintf('Creating data matrix...\n')
for i = 1:n_CAMELS_GB
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,n_CAMELS_GB)
    end
    
    t = datetime(CAMELS_GB_data.Q{i}(:,1),'ConvertFrom','datenum');
    Q = CAMELS_GB_data.Q{i}(:,2);
    P = CAMELS_GB_data.P{i}(:,2);
    PET = CAMELS_GB_data.PET{i}(:,2);
    
    t_mat{i} = t;
    Q_mat{i} = Q;
    P_mat{i} = P;
    PET_mat{i} = PET;
    
end

% We can now use the calculation function calc_Addor.m to calculate all the
% Addor et al. (2018) signatures, which are also contained in Coxon et al.
% (2020).
fprintf('Calculating signatures...\n')
CAMELS_GB_signatures = calc_Addor(Q_mat, t_mat, P_mat);
% Besides the signature values, the function also returns a list with 
% warnings and error messages. Most warnings come from our data check and 
% indicate that there are some NaN values in the time series. 
fprintf(CAMELS_GB_signatures.FDC_slope_error_str(1)+"\n")

% We can save the results as mat file which can be easily loaded into
% Matlab. Alternatively, we can save the results as txt file.
save(strcat(results_path,'CAMELS_GB_signatures.mat'),...
    '-struct','CAMELS_GB_signatures')
writetable(struct2table(CAMELS_GB_signatures),...
    strcat(results_path,'CAMELS_GB_signatures.txt'))

%% Compare TOSSH signatures to CAMELS signatures
% We can compare the signatures contained in CAMELS GB with the signatures
% calculated here to see if we get the same results.
makeScatterPlot(CAMELS_GB_signatures,CAMELS_GB_data,90)
saveFig(gcf,strcat('TOSSH_scatter_plot_GB'),fig_path,'-dpdf')
% Overall the results are very similar, but for some signatures there are
% large differences, which can be explained by different treatment of NaN
% values. The blue dots - which indicate a complete record - fall on a 1:1
% line, the other ones mostly don not. Note that some signatures are 
% extremly sensitive to NaN values, e.g. P-Q elasticity.

%% Calculation of new signatures
% We can also calculate some new signatures that are not provided with
% CAMELS GB.
% Note that some signatures take a while to calculate because there are
% many internal calculations (e.g. recession extration, fitting, etc.).

BaseflowRecessionK = NaN(n_CAMELS_GB,1);
BaseflowRecessionK_error_flag = NaN(n_CAMELS_GB,1);
BaseflowRecessionK_error_str = strings(n_CAMELS_GB,1);
RecessionParameters = NaN(n_CAMELS_GB,2);
RecessionParameters_error_flag = NaN(n_CAMELS_GB,1);
RecessionParameters_error_str = strings(n_CAMELS_GB,1);
EventRR = NaN(n_CAMELS_GB,1);
EventRR_error_flag= NaN(n_CAMELS_GB,1);
EventRR_error_str = strings(n_CAMELS_GB,1);

fprintf('Calculating new signatures...\n')
for i = 1:n_CAMELS_GB
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,n_CAMELS_GB)
    end
    
    % Since there are many missing values before 1989, we will only
    % consider the period from October 1989 to September 2009.
    t = datetime(CAMELS_GB_data.Q{i}(:,1),'ConvertFrom','datenum');
    Q = CAMELS_GB_data.Q{i}(:,2);
    P = CAMELS_GB_data.P{i}(:,2);
    PET = CAMELS_GB_data.PET{i}(:,2);
    
    indices = 1:length(t);
    start_ind = indices(t==datetime(1989,10,1));
    % in case time series starts after 1 Oct 1989
    if isempty(start_ind); start_ind = 1; end
    end_ind = indices(t==datetime(2009,9,30));
    t = t(start_ind:end_ind);
    Q = Q(start_ind:end_ind);
    P = P(start_ind:end_ind);
    PET = PET(start_ind:end_ind);    

    [BaseflowRecessionK(i),BaseflowRecessionK_error_flag(i),BaseflowRecessionK_error_str(i)] = ...
        sig_BaseflowRecessionK(Q,t,'recession_length',5);
    
    [RecessionParameters(i,:),~,RecessionParameters_error_flag(i),RecessionParameters_error_str(i)] ...
        = sig_RecessionAnalysis(Q,t,'fit_individual',false);
    
    [EventRR(i),EventRR_error_flag(i),EventRR_error_str(i)] = sig_EventRR(Q,t,P);
    
end
 
%% Plot maps
% We can plot the resulting signatures on a map using a plotting function
% Note that the mapping toolbox is needed for plotting.

plotMapUK(CAMELS_GB_data.gauge_lat,CAMELS_GB_data.gauge_lon,1./BaseflowRecessionK,...
    'attribute_name','K [d]','ID',CAMELS_GB_data.gauge_id,...
    'colour_scheme','bone','flip_colour_scheme',true,...
    'c_limits',[0 30],...
    'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','(a)','figure_name','RecessionK',...
    'save_figure',true,'figure_path',fig_path,'figure_type','-dpng')

plotMapUK(CAMELS_GB_data.gauge_lat,CAMELS_GB_data.gauge_lon,RecessionParameters(:,2),...
    'attribute_name','\beta [-]','ID',CAMELS_GB_data.gauge_id,...
    'colour_scheme','parula','flip_colour_scheme',true,...
    'c_limits',[1 2],...
    'c_lower_limit_open',true,'c_upper_limit_open',true,...
    'figure_title','(b)','figure_name','Recession_exponent',...
    'save_figure',true,'figure_path',fig_path,'figure_type','-dpng')

plotMapUK(CAMELS_GB_data.gauge_lat,CAMELS_GB_data.gauge_lon,EventRR,...
    'attribute_name','EventRR [-]','ID',CAMELS_GB_data.gauge_id,...
    'colour_scheme','pink','flip_colour_scheme',true,...
    'c_limits',[0.0 0.6],...
    'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','(c)','figure_name','EventRR',...
    'save_figure',true,'figure_path',fig_path,'figure_type','-dpng')

% We can see a clear patterns for all three signatures. The event runoff 
% mostly ratio follows climate arditiy. It is high along the very wet west
% coast and low in the drier south east. The recession constant and the
% recession exponent vary less smoothly in space as they are influenced 
% more strongly by the underlying geology. In areas underlain by highly 
% productive aquifers (e.g. Chalk), both the recession constant and the 
% recession exponent tend to be low.

%% Further analysis of resulting signatures
% We can also show more directly how certain signatures relate to certain
% catchment attributes by plotting them against each other. We use the
% CAMELS GB attribute "frac_high_perc" which quantifies the fraction of a
% catchment underlain by highly productive fractured aquifers to show the 
% influence of geology and aridity to show the influence of climate.

figure('pos',[100 100 400 300])
scatter(CAMELS_GB_data.aridity,BaseflowRecessionK,25,CAMELS_GB_data.frac_high_perc,'filled')
caxis([0 1]); c = colorbar; title(c,'% fract. aquifer [-]')
xlabel('Aridity [-]')
ylabel('K [1/d]')

figure('pos',[100 100 400 300])
scatter(CAMELS_GB_data.aridity,RecessionParameters(:,2),25,CAMELS_GB_data.frac_high_perc,'filled')
caxis([0 1]); c = colorbar; title(c,'% fract. aquifer [-]')
xlabel('Aridity [-]')
ylabel('\beta [-]'); ylim([0 3])

figure('pos',[100 100 400 300])
scatter(CAMELS_GB_data.aridity,EventRR,25,CAMELS_GB_data.frac_high_perc,'filled')
caxis([0 1]); c = colorbar; title(c,'% fract. aquifer [-]')
xlabel('Aridity [-]')
ylabel('EventRR [-]'); ylim([0 1])

% These plots correspond well with the patterns visible on the maps. The 
% event runoff ratio is strongly correlated with aridity. Geology 
% influences both the recession constant and the recession exponent, but 
% also the event runoff ratio, which is lower if a catchment is underlain 
% by a productive aquifer.

%% Further information
% Further information can be found in the online documentation:
% https://TOSSHtoolbox.github.io/TOSSH/ and in the other example scripts.
