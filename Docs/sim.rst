How to Run Simulations
========================

After downloading this repository, you can follow the follow steps to perform the stimulations. The main codes are based on MATLAB and have been verified in MATLAB R2019a. You can download all source codes from `this repository <https://github.com/pikipity/Stimulus-stimulus_transfer_based_on_time-frequency-joint_representation_SSVEP_BCIs.git>`_.

1. Download the Benchmark Dataset and put all 35 subjects' data in the `data` folder.

   + You can directly download data from the webpage of `Benchmark Dataset <http://bci.med.tsinghua.edu.cn/download.html>`_.
   + You also can use :file:`DownloadData.sh` to download data automatically. For linux, you can directly run this script. For windows, you can use `Cygwin <https://www.cygwin.com/>`_ to run this script. Please note that this script requires `wget <https://www.gnu.org/software/wget/>`_ and `7z <http://p7zip.sourceforge.net/>`_. You may refer the following commonds to run this script:
     
     .. code-block:: console

         cd <this_repository_path>
         sh ./DownloadData.sh

2. The following matlab toolboxes are required:

   + v2.1 MATLAB version of `Toolbox for Adaptive Fourier Decomposition <https://github.com/pikipity/Toolbox-for-Adaptive-Fourier-Decomposition>`_: The proposed stimulus-stimulus method requires this toolbox. If the version is lower than v2.1, the MAFD is not supported. So the version **MUST** be higher than or equal to v2.1.
   + To plot results, the following toolboxes are required. If you do not need to plot results, the following toolboxes are not needed.

     + `shadedErrorBar <https://github.com/raacampbell/shadedErrorBar>`_
     + `sigstar <https://github.com/raacampbell/sigstar>`_

3. In MATLAB, run :file:`Main.m` in the path of this repository. If the data and results are stored in other dictionary, you can change ``data_dir`` and ``figdir`` in :file:`Main.m`. 

.. note::

   The entire simulations may take a very **long** time and **large** RAM spaces. So you may need to run it step by step. 