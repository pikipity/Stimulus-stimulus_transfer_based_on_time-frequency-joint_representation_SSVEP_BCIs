Brief Introduction
=============================================================

Since the different stimuli are same in one period, we assume that the responses of our brains are also same. The following figure shows the fundamental assumption. The left polar plot shows the dynamic patterns of the signal states. For all stimuli, the patterns of the signal states are same. For different stimuli, the roration speeds of the signal states are different. Disturbed by the non-SSVEP-related components, the real EEG signals are generated as shown in the right side.

.. image:: https://github.com/pikipity/Stimulus-stimulus_transfer_based_on_time-frequency-joint_representation_SSVEP_BCIs/blob/main/images/plot_rotation.png?raw=true

From the transfer learning point of view, the phase domain can be considered as a high-level common domain that links sub-domains of different stimuli. The proposed method can be regarded as a asymmetric transformation. When the calibration data in the target domain is unavailable or not enough, the proposed method can be applied to generate the classification model of the target domain by transferring classification model of the non-target domain through the phase domain, which is illustrated in the following figure.

.. image:: https://github.com/pikipity/Stimulus-stimulus_transfer_based_on_time-frequency-joint_representation_SSVEP_BCIs/blob/main/images/TrnasferLearningDiagram.png?raw=true

In the proposed stimulus-stimulus transfer, the DP-MAFD is firstly applied to decompose averaged signals of non-target stimuli. Then, by analyzing decomposition coefficients, common coefficients can be found. Finally, these common coefficients can be adopted to construct SSVEP templates for target stimuli. The process is shown in the following figure.

.. image:: https://github.com/pikipity/Stimulus-stimulus_transfer_based_on_time-frequency-joint_representation_SSVEP_BCIs/blob/main/images/Plot_Progress_big.png?raw=true