.. _p4_usage_notes:

Usage Notes
========

Understanding Signature Parameters
----------------

Example parameters
^^^^^^^^^
Hydrologic signatures often need customized input parameters because rainfall patterns vary by region and time resolutions vary by dataset [1]_. 
The following lists example parameter files to help you get started. 
These parameters have been manually tuned for select signatures and U.S. gages at a **daily timestep**, using the CAMELS and HYSETS datasets in Caravan [2]_. 


1. Overlandflow signatures

- Parameter file: `params_overlandflow.csv <https://github.com/RY4GIT/TOSSH/tree/master/docs/_static/usages/params_overlandflow.csv>`_

- Use these parameters for the event separation algorithm `util_EventSeparation.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/utility_functions/util_EventSeparation.html>`_

- Parameters included: ``min_termination``, ``min_duration``, ``min_intensity_day``, ``min_intensity_day_during``, ``max_recessiondays`` per CONUS region

- In the table, ``ws_code`` column categorizes the CONUS regions using the first two digits of the USGS `gauge_id` from the CAMELS/HYSETS datasets. This classification is based on the downstream order system (see https://help.waterdata.usgs.gov/faq/sites/do-station-numbers-have-any-particular-meaning)

2. Recession-related signatures

- Parameter file: `params_recession.csv <https://github.com/RY4GIT/TOSSH/tree/master/docs/_static/usages/params_recession.csv>`_ 

- Use these parameters for the recession delineation algorithm `util_RecessionSegments.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/utility_functions/util_RecessionSegments.html>`_ 

- Parameters included: ``recession_length``, ``n_start``, ``eps``, ``filter_par`` for low or high rainfall area

**Note:** For tuning, we evaluated the parameter performance using approximately ten watersheds that have drainage area close to the median of each CONUS region. While these settings are a good starting point, please review and adjust them as needed for your specific application.

How to Adjust Parameters Yourself
^^^^^^^^^
If the provided parameters don't suit your dataset, here are some guidelines for tuning key parameters. 
We highly recommend using TOSSH visualization tools (`sig_EventGraphThresholds.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html>`_ and `sig_RecessionAnalysis.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_RecessionAnalysis.html>`_) to inspect how well the parameters fit your time series.

- Key parameters for `util_EventSeparation.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/utility_functions/util_EventSeparation.html>`_: 
	- ``min_termination``: 
		- This parameter controls how long the shortest storm event should be. The default value can be set to **48 hours** for daily data. With this, once-in-3-days rainfall is being separated as the new event. 
		- For persistent rainfall or slower-responding watersheds, increase the values (e.g., **72 hours**).
    - ``min_intensity_day`` and ``min_intensity_day_during``:
		- This parameter controls the minimum duration of the shortest storm event. The default value can be set to  **4.8 mm/day** for daily data (derived from 0.2 mm/hr × 24 hours). 
		- Adjust this value based on the region's rainfall intensity (e.g., use higher values for areas with significant rainfall). Selecting a small value for regions with high rainfall may result in small, noisy rainfall being detected at the start of an event.
    - ``max_recessiondays``:
        - The default value can be set to **8 days**  (no tuning is generally required).
        - Increasing this value helps capture long recession tails when no subsequent rainfall occurs. Allowing a larger ``max_recessiondays`` does no harm, as the ``min_termination`` parameter will start a new event when subsequent rainfall occurs.
    - ``min_duration``:
        - The default value can be set to **1 day** (no tuning is required for daily data).

- Key parameters for `util_RecessionSegments.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/utility_functions/util_RecessionSegments.html>`_:
    - ``recession_length``:
        - The default value can be set to **5 days** with ``eps=0.08`` for most cases.
        - For low-flow conditions, use longer thresholds (e.g., **10 days** with ``eps=0.01``).
        - For quick recessions after rainfall, reduce the value to **3 days** (but avoid including noisy recessions).

    - ``filter_par``:
        - The default value can be set to **0.925** (no tuning needed; used for smoothing).

    - ``n_start``:
        - The default value can be set to **0** for daily data. For hourly data, increase this value to mask flow immediately after the peak.

References
----------

.. [1] McMillan, H. K., Coxon, G., Araki, R., Salwey, S., Kelleher, C., Zheng, Y., et al. (2023). When good signatures go bad: Applying hydrologic signatures in large sample studies. Hydrological Processes, 37(9). https://doi.org/10.1002/hyp.14987

.. [2] Kratzert, F., Nearing, G., Addor, N., Erickson, T., Gauch, M., Gilon, O., Gudmundsson, L., Hassidim, A., Klotz, D., Nevo, S., Shalev, G., & Matias, Y. (2024). Caravan - A global community dataset for large-sample hydrology (1.4) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10968468