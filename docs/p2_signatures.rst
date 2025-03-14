.. _p2_signatures:

Signatures
==========

Below we provide a list of signatures with a short description organised according to the signature sets. 
A more detailed description of each signature can be found in the Matlab functions which can be opened by clicking on the function name.
At the end of the page we also provide a list with all remaining signatures contained in TOSSH which are not (yet) part of any set.

`Click here to browse through all Matlab functions contained in TOSSH. <./_static/matlab/TOSSH_code/index.html>`_

List of signature sets
----------------------

 
Basic set
^^^^^^^^^
The basic set of signatures (`calc_BasicSet.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_BasicSet.html>`_) 
is designed to cover the five components of a natural streamflow regime 
as defined by Richter et al. (1996) [1]_ and Poff et al. (1997) [2]_: 
magnitude, frequency, duration, timing, and rate of change.
As Poff et al. (1997) [2]_ state, these components “can be used to characterize the entire range of flows 
and specific hydrologic phenomena, such as floods or low flows, that are critical to the integrity of river ecosystems”. 

.. csv-table:: Basic set
   :file: ./_static/basic_set.csv
   :header-rows: 1
   :delim: ;

:sup:`a` Note that this signature function returns a value for each month.


Benchmark sets 
^^^^^^^^^^^^^^

Addor et al. (2018) set
"""""""""""""""""""""""
Addor et al. (2018) [3]_ use 15 signatures (`calc_Addor.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_Addor.html>`_) 
that “characterize different parts of the hydrograph, 
and [...] are sensitive to processes occurring over different time scales”. 
The signatures are selected from those commonly used in the literature, 
and are used to explore the strength of relationships between signatures and catchment attributes.

.. csv-table:: Addor et al. (2018)
   :file: ./_static/addor.csv
   :header-rows: 1
   :delim: ;


Sawicz et al. (2011) set
""""""""""""""""""""""""
Sawicz et al. (2011) [4]_ use 6 signatures (`calc_Sawicz.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_Sawicz.html>`_) 
drawn largely from Yadav et al. (2007) [5]_, 
that are chosen to be uncorrelated and to be linked to catchment function. 
The signatures are used to analyse hydrological similarity between catchments, 
and link the resulting clusters to climate and landscape attributes.

.. csv-table:: Sawicz et al. (2011)
   :file: ./_static/sawicz.csv
   :header-rows: 1
   :delim: ;


Euser et al. (2013) set 
"""""""""""""""""""""""
Euser et al. (2013) [6]_ use 8 signatures (`calc_Euser.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_Euser.html>`_) 
that represent different aspects of hydrological behaviour. 
The signatures are used to test the consistency of model performance, 
within the FARM model evaluation framework.
Some signatures are the same but applied to different parts of the time series, e.g. the low flow period (May to September) or the high flow period (November to April).

.. csv-table:: Euser et al. (2011)
   :file: ./_static/euser.csv
   :header-rows: 1
   :delim: ;

:sup:`a` Note that this signature function returns the whole flow duration curve.


Process-based signature sets 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

McMillan (2020) groundwater set
"""""""""""""""""""""""""""""""
Calculates 15 signatures (`calc_McMillan_Groundwater.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_McMillan_Groundwater.html>`_)
from McMillan (2020) [7]_, related to groundwater 
storage, groundwater dynamics and baseflow. These signatures come from 
previous experimental studies that link catchment or hillslope 
processes to streamflow response dynamics. Some signatures are 
implemented direct from the original papers, others are interpreted
from a qualitative description in the paper.

.. csv-table:: McMillan (2020) groundwater set
   :file: ./_static/groundwater.csv
   :header-rows: 1
   :delim: ;
   
:sup:`a` Note that these signature functions return multiple signatures.

:sup:`b` Note that this signature is a combination of two other signatures.


McMillan (2020) overland flow set
"""""""""""""""""""""""""""""""""
Calculates 10 overland flow (infiltration excess and saturation excess)
signatures (`calc_McMillan_OverlandFlow.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_McMillan_OverlandFlow.html>`_)
from McMillan (2020) [7]_. These signatures come from previous
experimental studies that link catchment or hillslope processes to
streamflow response dynamics. Some signatures are implemented direct
from the original papers, others are interpreted from a qualitative
description in the paper. 
Note that these signatures are all calculated with the same function.

.. csv-table:: McMillan (2020) overland flow set
   :file: ./_static/overland_flow.csv
   :header-rows: 1
   :delim: ;

Wu et al., (2021) overland flow set
"""""""""""""""""""""""""""""""""
A function `sig_EventGraphThresholds.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/signature_functions/sig_EventGraphThresholds.html>`_ contains 2 overland flow 
signatures from Wu et al., (2021) [8]_. The implementation in this toolbox is the adapted and coded version by Bolotin and McMillan (2025) [9]_. 

.. csv-table:: Wu et al., (2021) overland flow set
   :file: ./_static/overland_flow_Wu.csv
   :header-rows: 1
   :delim: ;


TOSSH: all remaining signatures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a list of all remaining signatures currently contained in TOSSH.

..
	(`calc_All.m <./_static/matlab/TOSSH_code/TOSSH/TOSSH_code/calculation_functions/calc_All.html>`_)


.. csv-table:: All remaining signatures
   :file: ./_static/remaining_signatures.csv
   :header-rows: 1
   :delim: ;
   
:sup:`a` Note that these signature functions return multiple signatures.

References
----------

.. [1] Richter, B.D., Baumgartner, J.V., Powell, J. and Braun, D.P., 1996. A method for assessing hydrologic alteration within ecosystems. Conservation biology, 10(4), pp.1163-1174.

.. [2] Poff, N.L., Allan, J.D., Bain, M.B., Karr, J.R., Prestegaard, K.L., Richter, B.D., Sparks, R.E. and Stromberg, J.C., 1997. The natural flow regime. BioScience, 47(11), pp.769-784.

.. [3] Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and Clark, M.P., 2018. A ranking of hydrological signatures based on their predictability in space. Water Resources Research, 54(11), pp.8792-8812.

.. [4] Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G., 2011. Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA. Hydrology and Earth System Sciences, 15(9), pp.2895-2911.

.. [5] Yadav, M., Wagener, T. and Gupta, H., 2007. Regionalization of constraints on expected catchment response behavior for improved predictions in ungauged basins. Advances in Water Resources, 30(8), pp.1756-1774.

.. [6] Euser, T., Winsemius, H.C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S. and Savenije, H.H.G., 2013. A framework to assess the realism of model structures using hydrological signatures. Hydrology and Earth System Sciences, 17(5), pp.1893-1912.

.. [7] McMillan, H.K., 2020. Linking hydrologic signatures to hydrologic processes: A review. Hydrological Processes, 34(6), pp.1393-1409.

.. [8] Wu, S., Zhao, J., Wang, H., & Sivapalan, M. (2021). Regional patterns and physical controls of streamflow generation across the conterminous United States. Water Resources Research, 57(6), e2020WR028086.

.. [9] Bolotin, L. A., & McMillan, H. (2024). A hydrologic signature approach to analysing wildfire impacts on overland flow. Hydrological Processes, 38(6). 

