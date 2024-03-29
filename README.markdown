# Stimulus-stimulus transfer based on time-frequency-joint representation in SSVEP-based BCI

This study is intended to achieve a general **stimulus-stimulus transfer** of SSVEP signals.

To achieve the stimulus-stimulus transfer, this study proposes a new perspective to analyze the SSVEP signals, which is the **time-frequency-joint representation**. In this representation, the SSVEP signals of different stimuli are synchronized according to the stimulus frequencies and phases. Then, the common components of different stimuli can be emphasied and become possible to be characterized and extracted. From the time-frequency-joint-representation point of view, a new SSVEP signal model is proposed for the stimulus-stimuls transfer. Moreover, an adaptive decomposition approach based on the multi-channel adaptive Fourier decomposition (MAFD) is designed to estimate components in the porposed SSVEP model, which is the MAFD with different phases (DP-MAFD). 

Codes in this repository implement the simulations for the Benchmark Dataset [1]. Simulation results show the classification performance of SSVEP templates constructed by the proposed stimulus-stimulus trnasfer. 

The related paper: Ze Wang et al., "Stimulus-stimulus transfer based on time-frequency-joint representation in SSVEP-based BCIs," *IEEE Trans. Biomed. Eng.*, 2022. DOI: [10.1109/TBME.2022.3198639](https://doi.org/10.1109/TBME.2022.3198639).

This repository follows the license [Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en).

## Citation

If you use this repository, please cite the related paper:

```
@article{wang_stimulus-stimulus_2022,
	title = {Stimulus-stimulus transfer based on time-frequency-joint representation in {SSVEP}-based {BCIs}},
	volume = {70},
	doi = {10.1109/TBME.2022.3198639},
	number = {2},
	journal = {IEEE Trans. Biomed. Eng.},
	author = {Wang, Ze and Wong, Chi Man and Rosa, Agostinho and Qian, Tao and Jung, Tzyy-Ping and Wan, Feng},
	year = {2022},
	pages = {603--615},
}
```

## Brief introduction to proposed stimulus-stimulus transfer

Since the different stimuli are same in one period, we assume that the responses of our brains are also same. The following figure shows the fundamental assumption. The left polar plot shows the dynamic patterns of the signal states. For all stimuli, the patterns of the signal states are same. For different stimuli, the roration speeds of the signal states are different. Disturbed by the non-SSVEP-related components, the real EEG signals are generated as shown in the right side.

![Example to show the fundamental assumption of the proposed stimulus-stimulus transfer](./images/plot_rotation.png)

From the transfer learning point of view, the phase domain can be considered as a high-level common domain that links sub-domains of different stimuli. The proposed method can be regarded as a asymmetric transformation. When the calibration data in the target domain is unavailable or not enough, the proposed method can be applied to generate the classification model of the target domain by transferring classification model of the non-target domain through the phase domain, which is illustrated in the following figure.

![Diagrams of the classification based on the proposed simulus-stimulus transfer without calibration data of target stimuli](./images/TrnasferLearningDiagram.png)

In the proposed stimulus-stimulus transfer, the DP-MAFD is firstly applied to decompose averaged signals of non-target stimuli. Then, by analyzing decomposition coefficients, common coefficients can be found. Finally, these common coefficients can be adopted to construct SSVEP templates for target stimuli. The process is shown in the following figure.

![Proposed stimulus-stimulus transfer](./images/Plot_Progress_big.png)

## Run simulations

After downloading this repository, you can follow the follow steps to perform the stimulations. The main codes are based on MATLAB and have been verified in MATLAB R2019a.

1. Download the Benchmark Dataset and put all 35 subjects' data in the `data` folder.

    + You can directly download data from the webpage of [Benchmark Dataset](http://bci.med.tsinghua.edu.cn/download.html).
    + You also can use `DownloadData.sh` to download data automatically. For linux, you can directly run this script. For windows, you can use [Cygwin](https://www.cygwin.com/) to run this script. Please note that this script requires [`wget`](https://www.gnu.org/software/wget/) and [`7z`](http://p7zip.sourceforge.net/). You may refer the following commonds to run this script. 

        ```
        cd <this_repository_path>
        sh ./DownloadData.sh
        ```

2. The following matlab toolboxes are required:

    + v2.1 MATLAB version of [Toolbox for Adaptive Fourier Decomposition](https://github.com/pikipity/Toolbox-for-Adaptive-Fourier-Decomposition): The proposed stimulus-stimulus method requires this toolbox. If the version is lower than v2.1, the MAFD is not supported. So the version **MUST** be higher than or equal to v2.1.
    + To plot results, the following toolboxes are required. If you do not need to plot results, the following toolboxes are not needed.
      + [shadedErrorBar](https://github.com/raacampbell/shadedErrorBar)
      + [sigstar](https://github.com/raacampbell/sigstar)

3. In MATLAB, run `Main.m` in the path of this repository. Note:
   + The entire simulations may take a very long time and a large RAM space. So you may want to run it step by step. 
   + If the data and results are stored in other dictionary, you can change `data_dir` and `figdir` in `Main.m`.

## Simulation results

In this simulation, the SSVEP templates generated by the proposed stimulus-stimulus method (without using calibration data of target stimuli) is compared with

+ SSVEP templates without using calibration data of target stimuli
  + Sine-cosine signals (CCA [2])
  + SSVEP templates generated by the superposition theory based SSVEP signal model (tlCCA [3])
+ SSVEP templates using calibration data of target stimuli
  + Averaged signals (eCCA [4])

### Comparisons of classification accuracy and ITR of porposed method and other methods without using calibration data of target stimuli

![Comparisons of ITR of porposed method and other methods without using calibration data of target stimuli](./ITR_Acc_Summary/ITR_summary_nocalibration_random_8.png)

![ITR_Comparisons of classification accuracy of porposed method and other methods without using calibration data of target stimuli](./ITR_Acc_Summary/Acc_summary_nocalibration_random_8.png)

### Comparisons of classification accuracy and ITR of porposed method and method using calibration data of target stimuli

![Comparisons of ITR of porposed method and method using calibration data of target stimuli](./ITR_Acc_Summary/ITR_summary_calibration_random_8.png)

![Comparisons of classification accuracy of porposed method and method using calibration data of target stimuli](./ITR_Acc_Summary/Acc_summary_calibration_random_8.png)

### Individual maximum ITR

![individual_ITR](./ITR_Acc_Summary/ITR_max_compare_all_8.png)

## References

1. Y. Wang, X. Chen, X. Gao, and S. Gao, "A benchmark dataset for SSVEP-based braincomputer interfaces," *IEEE Trans. Neural Syst. Rehabil. Eng.*, vol. 25, no. 10, pp. 1746–1752, 2017, doi: [10.1109/TNSRE.2016.2627556](https://doi.org/10.1109/TNSRE.2016.2627556).
2. Z. Lin, C. Zhang, W. Wu, and X. Gao, "Frequency recognition based on canonical correlation analysis for SSVEP-based BCIs," *IEEE Trans. Biomed. Eng.*, vol. 53, no. 12, pp. 2610–2614, 2006, doi: [10.1109/TBME.2006.886577](https://doi.org/10.1109/TBME.2006.886577).
3. C. M. Wong, Z. Wang, A. C. Rosa, C. L. P. Chen, T.-P. Jung, Y. Hu, and F. Wan, "Transferring subject-specific knowledge across stimulus
frequencies in SSVEP-based BCIs," *IEEE Trans. Automat. Sci. Eng.*, vol. 18, no. 2, pp. 552–563, 2021, doi: [10.1109/TASE.2021.3054741](https://doi.org/10.1109/TASE.2021.3054741).
4. M. Nakanishi, Y. Wang, Y.-T. Wang, Y. Mitsukura, and T.-P. Jung, "A high-speed brain speller using steady-state visual evoked potentials," *Int.
J. Neur. Syst.*, vol. 24, no. 06, p. 1450019, 2014, doi: [10.1142/S0129065714500191](https://doi.org/10.1142/S0129065714500191).


