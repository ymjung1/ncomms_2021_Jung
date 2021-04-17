# ncomms_2021_Jung

Codes are associated with the manuscript, Jung et al. "CD45 pre-exclusion from the tips of T cell microvilli prior to antigen recognition"

1. System requirements

MATLAB 
The provided scripts use Matlab (http://www.mathworks.com), version : '9.4.0.813654 (R2018a). No additional toolboxes are required.

2. Instructions for using the scripts with Demo

2.1. Expansion Airyscan analysis

2.1.1. Instruction
Run the 'run_airy_analysis.m' and select the 'Airy_zstack_example_5frames.mat' file  
2.1.2. Expected outcome
-Images display the ch1, ch2, and merged image, segmentation map, GP image, Histogram of GP distribution within the MV and CB areas, intensity, GP, SR plot for each frame. The result images are saved as a single 'avi' file. 
-An image displays plots of intensity for ch1 and ch2, GP, and SR as a function of distance from tips of MV calculated from the entire dataset.
2.1.3. Expected run time for the demo
80 sec 

2.2. STORM image analysis

2.2.1. Instruction
Run the 'run_storm_ONI.m'. Select the 'storm_example_ch0_Lsel.csv' file for ch0, and select the 'storm_example_ch1_CD45.csv' file for ch1. Modify threshold-parameters in the dialog box, for example, input '130,210' for 'X range', '180,260' for 'Y range' in the dialog box for setting a new ROI. 
2.2.2. Expected outcome
Data of the single molecules within the selected ROI for each channel, channel correction parameter, super-resoultion images for each channel, and the merged image between the two channels, Tip-locations determined by the segmentaion presented on a merged image, result of the 3D-Pair-correlaton (g(r)). 

2.2.3. Expected run time for the demo
70sec

2.3. Colocalization analysis
2.3.1. Instruction
Use the function 'colocalization_analysis.m.'.Loading the example data, 'colocalization_example_ch1.mat' for ch1 and 'colocalization_example_ch2.mat'
2.3.2. Expected outcome
Calculate the Pearson's correlation coefficient (R), Manders's overlap coefficient (MOC), Manders's correlation coefficients 1 and 2 (MCC1 and MCC2), coefficient kl and k2  [Reference: Manders, E. M. M., Verbeek, F. J. & Aten, J. A. Journal of Microscopy 169, 375-382, doi:10.1111/j.1365-2818.1993.tb03313.x (1993).]
2.3.1. Expected run time for the demo
0.2 sec 

3. Terms of use
Copyright (c) 2021 Yunmin Jung
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


4. Contact
Please contact yjung@lji.org for queries relating to the scripts.
