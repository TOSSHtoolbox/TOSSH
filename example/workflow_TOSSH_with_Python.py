### Example script showing how to run TOSSH from Python ###
# This script demonstrates how to run the TOSSH toolbox from Python using the MATLAB Engine API. It is based on the
# basic workflow (workflow 1) example provided in the TOSSH documentation.

### TOSSH workflow 1 - basic workflow (also shown in online documentation)
#
# This script shows the basic functionalities of TOSSH with some example data. The example data used in this workflow
# are taken from CAMELS-GB (Coxon et al., 2020), see README_example_data.txt for more information on data sources.
#
# Copyright (C) 2025
# This software is distributed under the GNU Public License Version 3.
# See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

# Check if required MATLAB toolboxes are installed. Required are:
#   - MATLAB (TOSSH was developed using Matlab R2020a)
#   - Statistics and Machine Learning Toolbox
#   - Optimization Toolbox

# The following Python packages are required:
#   - numpy
#   - pandas
#   - matplotlib

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# In addition, the Matlab engine for Python is required. This can be installed by navigating to the Matlab installation
# (e.g. .../MATLAB/R2024b/extern/engines/python/) and executing the setup.py script (e.g. python setup.py install).
# Users may want to install the Matlab engine in a virtual environment to avoid conflicts with other packages.
# More information can be found here:
# https://mathworks.com/help/matlab/matlab_external/python-setup-script-to-install-matlab-engine-api.html

import matlab.engine

eng = matlab.engine.start_matlab()  # Start the MATLAB engine
from functions.helper_TOSSH_with_python import run_tossh_function  # Load helper function that runs TOSSH signatures

# The code assumes your current working directory is the TOSSH directory. This can, however, be specified by the user.
mydir = os.getcwd()  # Get current working directory, e.g. for loading data in subfolders
sys.path.append(mydir)  # Add current working directory to path (optional)

# Input data
# Every signature requires a streamflow (Q) time series and a corresponding date vector (t). Here we use Python's
# datetime module to create a date and covert it internally before passing it to Matlab (in the run_tossh_function).
# Here is an example of Q and t vectors in the correct format:
Q = np.array([1.14, 1.07, 2.39, 2.37, 1.59, 1.25, 1.35, 1.16, 1.27, 1.14])
t = pd.date_range(start='1999-10-01', periods=10, freq='D')

# Typically, users will have their own data which they want to analyse. We provide an example file to get a more
# realistic time series. The example file also contains precipitation (P), potential evapotranspiration (PET), and
# temperature (T) data, which are required for some signatures.
path = mydir + '/example_data/'  # specify path

# The original TOSSH examples use .mat files, but here we will load the data from a csv file (which has been added).
data = pd.read_csv(path + '33029_daily.csv', delimiter=',')
t = pd.to_datetime(data['t']).values
Q = data['Q'].values  # streamflow [mm/day]
P = data['P'].values  # precipitation [mm/day]
# Note: PET and T example data are provided but not used here.
# PET = data['PET'].values  # potential evapotranspiration [mm/day]
# T = data['T'].values  # temperature [degC]
# Note that the data format is assumed to be an array (for t a datetime64 array and for the other data a float array).
# Other formats may also be used (e.g. lists, pandas Series/DataFrames), but would require adjustments in the code.

# Plot data
# We can plot the data to get a first idea of the hydrograph.
plt.figure(figsize=(7, 4))
plt.plot(t, Q, 'k-', linewidth=1.0)
plt.xlabel('Date')
plt.gca().xaxis.set_major_locator(plt.MaxNLocator(5))
plt.ylabel('Streamflow [mm/day]')
plt.show()

# More information on the input data can be found here:
# https://TOSSHtoolbox.github.io/TOSSH/p1_overview.html.

# Calculate signatures
# Once the input data are loaded, we can calculate different signatures. We start by calculating the mean flow Q_mean.
# Note that the following function call only returns the Q_mean value and not other potential function outputs, because
# nargout is set to 1 by default. We will later show how more outputs can be retrieved.
Q_mean = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_Q_mean', eng, Q, t)
print(f'Q_mean = {Q_mean:.2f}')

# Some signatures can be calculated using different methods and/or parameter values. For example, there are different
# options to calculate the baseflow index (BFI). The default method is the UKIH smoothed minima method with a parameter
# of 5 days.
BFI_UKIH = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_BFI', eng, Q, t)
print(f'BFI_UKIH = {BFI_UKIH:.2f}')
# Alternatively, we can use the Lyne-Hollick filter with a filter parameter
# of 0.925.
BFI_LH = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_BFI', eng, Q, t, method='Lyne_Hollick')
print(f'BFI_LH = {BFI_LH:.2f}')
# We can also change the parameter value of the UKIH method to 10 days.
BFI_UKIH10 = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_BFI', eng, Q, t, method='UKIH', parameters=10)
print(f'BFI_UKIH10 = {BFI_UKIH10:.2f}')
# As we can see, all three options lead to slightly different values. More details and examples on the different
# methods/parameters can be found in the code of each function (e.g. sig_BFI.m).

# Some signatures also require precipitation (P) input time series. For example, the total runoff ratio requires both Q
# and P time series.
TotalRR = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_TotalRR', eng, Q, t, P)
print(f'TotalRR = {TotalRR:.2f}')

# Some signature functions come with a plotting functionality. For example, we can calculate the slope of the flow
# duration curve (FDC) and plot the result.
FDC_slope = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_FDC_slope', eng, Q, t, plot_results=True)
print(f'FDC_slope = {FDC_slope:.2f}')

# Some signatures are combinations of existing signatures, e.g. the baseflow fraction (K_b) defined as the ratio between
# mean baseflow Q_b and mean precipitation P. This signature can also be calculated as K_b = BFI*TotalRR. We therefore
# do not provide an extra signature function, but suggest to use the two existing functions.
K_b = (run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_BFI', eng, Q, t) *
       run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_TotalRR', eng, Q, t, P))
print(f'K_b = {K_b:.2f}')

# More information on the signatures contained in TOSSH can be found here:
# https://sebastiangnann.github.io/TOSSH_development/p2_signatures.html

# Warnings and errors
# Each signature function can return a warning/error output. These warning/error outputs indicate problems during
# signature calculation, but they do not stop code execution like a normal Matlab error would do. Two outputs can be
# retrieved: an error flag (error_flag) that corresponds to a certain type of warning/error, and a string (error_str)
# that decribes the warning/error. If multiple warnings/errors occur, they are all listed in the error string, starting
# with the one that occurred last.

# A warning (error_flag = 1) typically indicates that the signature can be calculated but should be interpreted with
# care, e.g. because there are NaN values in the time series. To retrieve the error flags and other outputs beyond the
# signature values, we need to specify the number of output arguments (nargout) in the run_tossh_function. For example,
# in the following it is set to 3, so that we can retrieve the error flag and error string. Functions that return even
# more outputs can be handled in the same way by setting nargout accordingly.
Q[:10] = np.nan
Q_mean, error_flag, error_str = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_Q_mean', eng, Q, t, nargout=3)
print(f'Q_mean = {Q_mean:.2f}')
print(error_str)
# We get the same mean value as before since the ten removed values do not influence the result much. In other cases,
# significant amounts of NaN entries might cause more problems.

# An error (error_flag > 1) indicates that the signature could not be calculated, e.g. because there is a problem with
# the input data. For example, if the input time series contains negative and thus physically impossible values, NaN is
# returned.
Q[:10] = -1.0
Q_mean, error_flag, error_str = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_Q_mean', eng, Q, t, nargout=3)
print(f'Q_mean = {Q_mean:.2f}')
print(error_str)

# Since these warnings/errors do not stop the execution of the code, we can run the signature code for many catchments
# without breaking, even if for some of the catchments the signature cannot be calculated.

# There are also "normal" errors which can happen if the input parameters are specified incorrectly (wrong format, wrong
# range, etc.). For example, if we swap Q and t in the input, we will get an error when directly working in Matlab. In
# this Python workflow, this is not the case due to the way we internally convert the input data (e.g. to datenum).
# However, the result will be wrong, as the function will interpret the input incorrectly.
Q_mean, error_flag, error_str = run_tossh_function(r'D:/Matlab/TOSSH/TOSSH_code', 'sig_Q_mean', eng, t, Q, nargout=3)
print(f'Q_mean = {Q_mean:.2f}')
print(error_str)
# Such errors should be easily avoidable by specifying all inputs correctly.

# Further information
# Further information can be found in the online documentation: https://TOSSHtoolbox.github.io/TOSSH/
