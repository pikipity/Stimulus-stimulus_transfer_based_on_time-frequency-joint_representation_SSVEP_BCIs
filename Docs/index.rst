.. Stimulus-stimulus_transfer_based_on_time-frequency-joint_representation_SSVEP_BCIs documentation master file, created by
   sphinx-quickstart on Thu Aug 11 08:59:02 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Stimulus-stimulus Transfer Based on Time-frequency-joint Representation in SSVEP-based BCIs
==============================================================================================================

.. toctree::
   :maxdepth: 1
   
   intro
   sim
   res
   license

This study is intended to achieve a general **stimulus-stimulus transfer** of SSVEP signals.

To achieve the stimulus-stimulus transfer, this study proposes a new perspective to analyze the SSVEP signals, which is the **time-frequency-joint representation**. In this representation, the SSVEP signals of different stimuli are synchronized according to the stimulus frequencies and phases. Then, the common components of different stimuli can be emphasied and become possible to be characterized and extracted. From the time-frequency-joint-representation point of view, a new SSVEP signal model is proposed for the stimulus-stimuls transfer. Moreover, an adaptive decomposition approach based on the multi-channel adaptive Fourier decomposition (MAFD) is designed to estimate components in the porposed SSVEP model, which is the MAFD with different phases (DP-MAFD). 

Codes in `this repository <https://github.com/pikipity/Stimulus-stimulus_transfer_based_on_time-frequency-joint_representation_SSVEP_BCIs.git>`_ implement the simulations for the Benchmark Dataset [1]. Simulation results show the classification performance of SSVEP templates constructed by the proposed stimulus-stimulus trnasfer. 

The related paper: Ze Wang et al., "Stimulus-stimulus transfer based on time-frequency-joint representation in SSVEP-based BCIs," *IEEE Trans. Biomed. Eng.*, 2022. DOI: `10.1109/TBME.2022.3198639 <https://doi.org/10.1109/TBME.2022.3198639>`_.

References
--------------

1. Y. Wang, X. Chen, X. Gao, and S. Gao, "A benchmark dataset for SSVEP-based braincomputer interfaces," *IEEE Trans. Neural Syst. Rehabil. Eng.*, vol. 25, no. 10, pp. 1746-1752, 2017, doi: `10.1109/TNSRE.2016.2627556 <https://doi.org/10.1109/TNSRE.2016.2627556>`_.



