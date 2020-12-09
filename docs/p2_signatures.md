---
title: Signatures
nav: true
---

- 
{:toc}

Below you can find a list of signatures with a short description organised according to the signature sets. 
A more detailed description of each signature can be found in the Matlab functions which can be opened by clicking on the function name.
At the end of the page you can also find a list with all remaining signatures contained in TOSSH which are not (yet) part of any set.

<a href="./matlab/TOSSH_code/index.html" target="_blank">Click here to browse through all Matlab functions contained in TOSSH.</a>
 

## List of signature sets

 
### Basic set 
The basic set of signatures (<a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_BasicSet.html" target="_blank">calc_BasicSet.m</a>) 
are designed to cover the five components of a natural streamflow regime 
as defined by Richter et al. (1996) and Poff et al. (1997): 
magnitude, frequency, duration, timing and rate of change.
As Poff et al. (1997) state, these components “can be used to characterize the entire range of flows 
and specific hydrologic phenomena, such as floods or low flows, that are critical to the integrity of river ecosystems”. 

| Signature | Function (click for code) | Unit | Description&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 
| ----------- | ----------- | ----------- | ----------- | 
| **Magnitude** |
| Q_mean | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Q_mean.html" target="_blank">sig_Q_mean.m</a> | mm/timestep | Mean streamflow |
| Q5 | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_percentile.html" target="_blank">sig_x_percentile.m</a> | mm/timestep | 5-th streamflow percentile |
| Q95 | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_percentile.html" target="_blank">sig_x_percentile.m</a> | mm/timestep | 95-th streamflow percentile |
| Q_mean_monthly | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Q_mean_monthly.html" target="_blank">sig_Q_mean_monthly.m</a> | mm/timestep | Mean monthly streamflow |
| Q_7_day_min | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Q_n_day_min.html" target="_blank">sig_Q_n_day_min.m</a> | mm/timestep | 7-day minimum streamflow |
| BFI | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_BFI.html" target="_blank">sig_BFI.m</a> | - | Baseflow index |
| CoV | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_CoV.html" target="_blank">sig_CoV.m</a> | - | Coefficient of variation |
| **Frequency** |
| high_Q_frequency | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_frequency.html" target="_blank">sig_x_Q_frequency.m</a> | - | High flow frequency |
| low_Q_frequency | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_frequency.html" target="_blank">sig_x_Q_frequency.m</a> | - | Low flow frequency |
| zero_Q_frequency | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_frequency.html" target="_blank">sig_x_Q_frequency.m</a> | - | Zero flow frequency |
| **Duration** |
| high_Q_duration | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_duration.html" target="_blank">sig_x_Q_duration.m</a> | timestep | High flow duration |
| low_Q_duration | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_duration.html" target="_blank">sig_x_Q_duration.m</a> | timestep | Low flow duration |
| zero_Q_duration | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_duration.html" target="_blank">sig_x_Q_duration.m</a> | timestep | Zero flow duration |
| **Timing** |
| HFD_mean | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_HFD_mean.html" target="_blank">sig_HFD_mean.m</a> | day of year | Half flow date |
| HFI_mean | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_HFI_mean.html" target="_blank">sig_HFI_mean.m</a> | days | Half flow interval |
| **Rate of change** |
| AC1 | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Autocorrelation.html" target="_blank">sig_Autocorrelation.m</a> | - | Lag-1 autocorrelation |
| FDC_slope | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_FDC_slope.html" target="_blank">sig_FDC_slope.m</a> | - | Slope of the flow duration curve |
| BaseflowRecessionK | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_BaseflowRecessionK.html" target="_blank">sig_BaseflowRecessionK.m</a> | 1/d | Recession constant |


### Benchmark sets 

#### Addor et al. (2018) set 
Addor et al. (2018) use 15 signatures (<a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_Addor.html" target="_blank">calc_Addor.m</a>)
that “characterize different parts of the hydrograph, 
and [...] are sensitive to processes occurring over different time scales”. 
The signatures were selected from those commonly used in the literature, 
and are used to explore the strength of relationships between signatures and catchment attributes.

| Signature | Function (click for code) | Unit | Description&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 
| ----------- | ----------- | ----------- | ----------- |  
| Q_mean | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Q_mean.html" target="_blank">sig_Q_mean.m</a> | mm/timestep | Mean streamflow |
| TotalRR | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_TotalRR.html" target="_blank">sig_TotalRR.m</a> | - | Total runoff ratio |
| QP_elasticity | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_QP_elasticity.html" target="_blank">sig_QP_elasticity.m</a> | - | Streamflow-precipitation elasticity |
| FDC_slope | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_FDC_slope.html" target="_blank">sig_FDC_slope.m</a> | - | Slope of the flow duration curve |
| BFI | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_BFI.html" target="_blank">sig_BFI.m</a> | - | Baseflow index |
| HFD_mean | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_HFD_mean.html" target="_blank">sig_HFD_mean.m</a> | day of year | Half flow date |
| Q5 | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_percentile.html" target="_blank">sig_x_percentile.m</a> | mm/timestep | 5th streamflow percentile (low flows) |
| Q95 | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_percentile.html" target="_blank">sig_x_percentile.m</a> | mm/timestep | 95th streamflow percentile (high flows) |
| high_Q_freq | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_frequency.html" target="_blank">sig_x_Q_frequency.m</a> | - | High flow frequency |
| high_Q_dur | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_duration.html" target="_blank">sig_x_Q_duration.m</a> | timestep | High flow duration |
| low_Q_freq | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_frequency.html" target="_blank">sig_x_Q_frequency.m</a> | - | Low flow frequency |
| low_Q_dur | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_duration.html" target="_blank">sig_x_Q_duration.m</a> | timestep | Low flow duration |
| zero_Q_freq | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_x_Q_frequency.html" target="_blank">sig_x_Q_frequency.m</a> | - | Zero flow frequency |


#### Sawicz et al. (2011) set 
Sawicz et al. (2011) use 6 signatures (<a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_Sawicz.html" target="_blank">calc_Sawicz.m</a>) 
drawn largely from Yadav et al. (2007), 
that are chosen to be uncorrelated and to be linked to catchment function. 
The signatures are used to analyse hydrological similarity between catchments, 
and link the resulting clusters to climate and landscape attributes.

| Signature | Function (click for code) | Unit | Description&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |  
| ----------- | ----------- | ----------- | ----------- | 
| TotalRR | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_TotalRR.html" target="_blank">sig_TotalRR.m</a> | - | Total runoff ratio |
| FDC_slope | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_FDC_slope.html" target="_blank">sig_FDC_slope.m</a> | - | Slope of the flow duration curve |
| BFI | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_BFI.html" target="_blank">sig_BFI.m</a> | - | Baseflow index |
| QP_elasticity | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_QP_elasticity.html" target="_blank">sig_QP_elasticity.m</a> | - | Streamflow-precipitation elasticity |
| SnowDayRatio | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_SnowDayRatio.html" target="_blank">sig_SnowDayRatio.m</a> | - | Snow day ratio |
| RLD | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_RisingLimbDensity.html" target="_blank">sig_RisingLimbDensity.m</a> | 1/timestep | Rising limb density |

#### Euser et al. (2013) set 
Euser et al. (2013) use 8 signatures (<a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_Euser.html" target="_blank">calc_Euser.m</a>) 
that represent different aspects of hydrological behaviour. 
The signatures are used to test the consistency of model performance, 
within the FARM model evaluation framework.

| Signature | Function (click for code) | Unit | Description&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 
| ----------- | ----------- | ----------- | ----------- | 
| AC1 | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Autocorrelation.html" target="_blank">sig_Autocorrelation.m</a> | - | Lag-1 autocorrelation |
| AC1_low | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Autocorrelation.html" target="_blank">sig_Autocorrelation.m</a> | - | Lag-1 autocorrelation for low flow period |
| RLD | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_RisingLimbDensity.html" target="_blank">sig_RisingLimbDensity.m</a> | 1/timestep | Rising limb density |
| PeakDistribution | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_PeakDistribution.html" target="_blank">sig_PeakDistribution.m</a> | - | Slope of distribution of peaks |
| PeakDistribution_low | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_PeakDistribution.html" target="_blank">sig_PeakDistribution.m</a> | - | Slope of distribution of peaks for low flow period |
| FDC | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_FDC.html" target="_blank">sig_FDC.m</a> | - | Flow duration curve |
| FDC_low | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_FDC.html" target="_blank">sig_FDC.m</a> | - | Flow duration curve for low flow period |
| FDC_high | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_FDC.html" target="_blank">sig_FDC.m</a> | - | Flow duration curve for high flow period |


###  Process-based signature sets 

#### McMillan (2020) groundwater set 
Calculates 15 signatures (<a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_McMillan_Groundwater.html" target="_blank">calc_McMillan_Groundwater.m</a>)
from McMillan (2020), related to groundwater 
storage, groundwater dynamics and baseflow. These signatures come from 
previous experimental studies that link watershed or hillslope 
processes to streamflow reponse dynamics. Some signatures are 
implemented direct from the original papers, others are interpreted
from a qualitative description in the paper.

| Signature | Function (click for code) | Unit | Description&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|
| ----------- | ----------- | ----------- | ----------- | 
| TotalRR | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_TotalRR.html" target="_blank">sig_TotalRR.m</a> | - | Total runoff ratio |
| RR_Seasonality | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_RR_Seasonality.html" target="_blank">sig_RR_Seasonality.m</a> | - | Runoff ratio seasonality |
| EventRR | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventRR.html" target="_blank">sig_EventRR.m</a> | - | Event runoff ratio |
| StorageFraction | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_StorageFraction.html" target="_blank">sig_StorageFraction.m</a> | - | Ratio between active and total storage |
| Recession_a_Seasonality | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_SeasonalVarRecessions.html" target="_blank">sig_SeasonalVarRecessions.m</a> | - | Seasonal variations in recession rate |
| AverageStorage | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_StorageFromBaseflow.html" target="_blank">sig_StorageFromBaseflow.m</a> | - | Average storage from average baseflow and storage-discharge relationship|
| RecessionParameters | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_RecessionAnalysis.html" target="_blank">sig_RecessionAnalysis.m</a> | - | Recession analysis parameters approximate storage-discharge relationship |
| MRC_num_segments | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_MRC_SlopeChanges.html" target="_blank">sig_MRC_SlopeChanges.m</a> | - | Number of different segments in MRC |
| First_Recession_Slope | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_MRC_SlopeChanges.html" target="_blank">sig_MRC_SlopeChanges.m</a> | mm/timestep<sup>2</sup> | Steep section of the master recession curve = storage that is quickly depleted |
| Mid_Recession_Slope | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_MRC_SlopeChanges.html" target="_blank">sig_MRC_SlopeChanges.m</a> | mm/timestep<sup>2</sup> | Mid section of the master recession curve = water retention capacity of the watershed |
| Spearmans_rho | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_RecessionUniqueness.html" target="_blank">sig_RecessionUniqueness.m</a> | - | Non-uniqueness in the recession relationship |
| EventRR_TotalRR_ratio<sup>1</sup> | EventRR/TotalRR | - | Ratio between event and total runoff ratio |
| VariabilityIndex | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_VariabilityIndex.html" target="_blank">sig_VariabilityIndex.m</a> | - | Variability index of flow |
| BFI | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_BFI.html" target="_blank">sig_BFI.m</a> | - | Baseflow index |
| BaseflowRecessionK | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_BaseflowRecessionK.html" target="_blank">sig_BaseflowRecessionK.m</a> | 1/d | Recession constant |

<sup>1</sup>Note that this signature is a combination of two other signatures.


#### McMillan (2020) overland flow set
Calculates 10 overland flow (infiltration excess and saturation excess)
signatures (<a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_McMillan_OverlandFlow.html" target="_blank">calc_McMillan_OverlandFlow.m</a>)
from McMillan (2020). These signatures come from previous
experimental studies that link watershed or hillslope processes to
streamflow reponse dynamics. Some signatures are implemented direct
from the original papers, others are interpreted from a qualitative
description in the paper. 

| Signature | Function (click for code) | Unit | Description&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 
| ----------- | ----------- | ----------- | ----------- | 
| IE_effect | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | - | Infiltration excess importance |
| SE_effect | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | - | Saturation excess importance |
| IE_thresh_signif | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | - | Infiltration excess threshold significance (in a plot of quickflow volume vs. maximum intensity) | 
| IE_thresh | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | mm/timestep | Infiltration excess threshold location (in a plot of quickflow volume vs. maximum intensity) | 
| SE_thresh_signif | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | - | Saturation excess threshold significance (in a plot of quickflow volume vs. total precipitation) |
| SE_thresh | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | mm | Saturation excess threshold location (in a plot of quickflow volume vs. total precipitation) |
| SE_slope | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | mm/mm | Saturation excess threshold above-threshold slope (in a plot of quickflow volume vs. total precipitation) |
| Storage_thresh_signif | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | - | Storage/saturation excess threshold significance (in a plot of quickflow volume vs. antecedent precipitation index + total precipitation) | 
| Storage_thresh | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | mm | storage/saturation excess threshold location (in a plot of quickflow volume vs. antecedent precipitation index + total precipitation) |
| min_Qf_perc | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html" target="_blank">sig_EventGraphThresholds.m</a> | % | Minimum quickflow as a percentage of precipitation | 


### TOSSH: all remaining signatures 

This is a list of all remaining signatures 
(<a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_All.html" target="_blank">calc_All.m</a>)
currently contained in TOSSH. 

| Signature | Function (click for code) | Unit | Description&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 
| ----------- | ----------- | ----------- | ----------- |
| BaseflowMagnitude | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_BaseflowMagnitude.html" target="_blank">sig_BaseflowMagnitude.m</a> | mm | Difference between maximum and minimum of annual baseflow regime |
| ResponseTime | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_ResponseTime.html" target="_blank">sig_ResponseTime.m</a> | timestep | Catchment response time |
| FlashinessIndex | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_FlashinessIndex.html" target="_blank">sig_FlashinessIndex.m</a> | - | Richards-Baker flashiness idex |
| HFI_mean | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_HFI_mean.html" target="_blank">sig_HFI_mean.m</a> | timestep | Half flow interval |
| PQ_Curve<sup>1</sup> | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_PQ_Curve.html" target="_blank">sig_PQ_Curve.m</a> | - | Slopes and breakpoints in cumulative P-Q regime curve |
| Q_n_day_max | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Q_n_day_max.html" target="_blank">sig_Q_n_day_max.m</a> | mm/timestep | n-day maximum streamflow |
| Q_skew | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Q_skew.html" target="_blank">sig_Q_skew.m</a> | mm<sup>3</sup>/timestep<sup>3</sup> | Skewness of streamflow|
| Q_var | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_Q_var.html" target="_blank">sig_Q_var.m</a> | mm<sup>2</sup>/timestep<sup>2</sup> | Variance of streamflow |
| RecessionK_part | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_RecessionParts.html" target="_blank">sig_RecessionParts.m</a> | 1/timestep | Recession constant of early/late (exponential) recessions |
| SeasonalTranslation<sup>1</sup>| <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_SeasonalTranslation.html" target="_blank">sig_SeasonalTranslation.m</a> | - | Amplitude ratio and phase shift between seasonal forcing and flow cycles  |
| SnowStorage | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_SnowStorage.html" target="_blank">sig_SnowStorage.m</a> | mm | Snow storage derived from cumulative P-Q regime curve |
| StorageFraction<sup>1</sup> | <a href="./matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_StorageFraction.html" target="_blank">sig_StorageFraction.m</a> | - | Ratio of active to total storage |

<sup>1</sup>Note that these signature functions return multiple signatures.

## References

Poff, N.L., Allan, J.D., Bain, M.B., Karr, J.R., Prestegaard, K.L., Richter, B.D., Sparks, R.E. and Stromberg, J.C., 1997. The natural flow regime. BioScience, 47(11), pp.769-784.

Richter, B.D., Baumgartner, J.V., Powell, J. and Braun, D.P., 1996. A method for assessing hydrologic alteration within ecosystems. Conservation biology, 10(4), pp.1163-1174.

Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and Clark, M.P., 2018. A ranking of hydrological signatures based on their predictability in space. Water Resources Research, 54(11), pp.8792-8812.

Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G., 2011. Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA. Hydrology & Earth System Sciences Discussions, 8(3).

Yadav, M., Wagener, T. and Gupta, H., 2007. Regionalization of constraints on expected watershed response behavior for improved predictions in ungauged basins. Advances in water resources, 30(8), pp.1756-1774.

Euser, T., Winsemius, H.C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S. and Savenije, H.H.G., 2013. A framework to assess the realism of model structures using hydrological signatures. Hydrology and Earth System Sciences, 17 (5), 2013.

McMillan, H., 2020. Linking hydrologic signatures to hydrologic processes: A review. Hydrological Processes, 34(6), pp.1393-1409.

