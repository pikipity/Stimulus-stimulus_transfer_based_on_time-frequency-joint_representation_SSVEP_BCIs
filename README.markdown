# Phase-domain transformation based stimulus-stimulus transfer for SSVEP-based BCI

This study is intended to achieve a general **stimulus-stimulus transfer** of SSVEP signals.

To achieve the stimulus-stimulus transfer, this study proposes a new perspective to analyze the SSVEP signals, which is the **phase domain**. In the phase domain, the SSVEP signals of different stimuli are synchronized according to the stimulus frequencies and phases. Then, the common components of different stimuli can be emphasied and become possible to be characterized and extracted. From the phase-domain point of view, a new SSVEP signal model is proposed for the stimulus-stimuls transfer. Moreover, an adaptive decomposition approach based on the multi-channel adaptive Fourier decomposition (MAFD) is designed to estimate components in the porposed SSVEP model, which is the MAFD with different phases (DP-MAFD). 

Codes in this repository implement the simulations for the Benchmark Dataset. Simulation results show the classification performance of SSVEP templates constructed by the proposed stimulus-stimulus trnasfer. 

## Brief introduction to proposed stimulus-stimulus transfer

Since the different stimuli are same in one period, we assume that the responses of our brains are also same. The following figure shows the fundamental assumption. The left polar plot shows the dynamic patterns of the signal states. For all stimuli, the patterns of the signal states are same. For different stimuli, the roration speeds of the signal states are different. Disturbed by the non-SSVEP-related components, the real EEG signals are generated as shown in the right side.

![Example to show the fundamental assumption of the proposed stimulus-stimulus transfer](./images/plot_rotation.png)

From the transfer learning point of view, the phase domain can be considered as a high-level common domain that links sub-domains of different stimuli. The proposed method can be regarded as a asymmetric transformation. When the calibration data in the target domain is unavailable or not enough, the proposed method can be applied to generate the classification model of the target domain by transferring classification model of the non-target domain through the phase domain, which is illustrated in the following figure.

![Diagrams of the classification based on the proposed simulus-stimulus transfer without calibration data of target stimuli](./images/TrnasferLearningDiagram.png)

In the proposed stimulus-stimulus transfer, the DP-MAFD is firstly applied to decompose averaged signals of non-target stimuli. Then, by analyzing decomposition coefficients, common coefficients can be found. Finally, these common coefficients can be adopted to construct SSVEP templates for target stimuli. The process is shown in the following figure.

![Proposed stimulus-stimulus transfer](./images/Plot_Progress_big.png)

## Codes of simulations

After downloading this repository, you can follow the follow steps to perform the stimulations. The main codes are based on MATLAB and have been verified in MATLAB R2019a

1. Download the Benchmark Dataset and put all 35 subjects' data in the `data` folder.

    + You can directly download data from the webpage of [Benchmark Dataset](http://bci.med.tsinghua.edu.cn/download.html).
    + You also can use `DownloadData.sh` to download data automatically. For linux, you can directly run this script. For windows, you can use [Cygwin](https://www.cygwin.com/) to run this script. Please note that this script requires [`wget`](https://cygwin.com/packages/summary/wget.html) and [`7z`](https://cygwin.com/packages/summary/p7zip.html). You may refer the following commonds to run this script. 

        ```
        cd <this_repository_path>
        sh ./DownloadData.sh
        ```

2. This simulation requires the following matlab toolboxes:

    + v2.1 MATLAB version of [Toolbox for Adaptive Fourier Decomposition](https://github.com/pikipity/Toolbox-for-Adaptive-Fourier-Decomposition): The proposed stimulus-stimulus method requires this toolbox. If the version is lower than v2.1, the MAFD is not supported. So the version must be higher than or equal to v2.1.
    + To plot results, the following toolboxes are required. If you do not need to plot results, the following toolboxes are not needed.
      + [shadedErrorBar](https://github.com/raacampbell/shadedErrorBar)
      + [sigstar](https://github.com/raacampbell/sigstar)

3. In MATLAB, run `Main.m` in the path of this repository. 

## Simulation results

### Comparisons of classification accuracy and ITR of porposed method and other methods without using calibration data of target stimuli

![Comparisons of ITR of porposed method and other methods without using calibration data of target stimuli](./ITR_Acc_Summary/ITR_summary_nocalibration_random_8.png)

![ITR_Comparisons of classification accuracy of porposed method and other methods without using calibration data of target stimuli](./ITR_Acc_Summary/Acc_summary_nocalibration_random_8.png)

### Comparisons of classification accuracy and ITR of porposed method and method using calibration data of target stimuli

![Comparisons of ITR of porposed method and method using calibration data of target stimuli](./ITR_Acc_Summary/ITR_summary_calibration_random_8.png)

![Comparisons of classification accuracy of porposed method and method using calibration data of target stimuli](./ITR_Acc_Summary/Acc_summary_calibration_random_8.png)

### Individual maximum ITR

![individual_ITR](./ITR_Acc_Summary/ITR_max_compare_all_8.png)

