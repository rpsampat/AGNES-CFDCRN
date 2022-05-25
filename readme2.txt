This folder contains the code used and adjusted by Maaike de Wit as part of her MSc thesis in Aerospace Engineering at TU Delft and builds upon the work performed by Rishikesh Sampat.

The code still has the same overall package set up as described in 'readme.txt' and important changes with respect to the original code were logged in 'CHANGELOG.txt'.

The thesis for which adjustments were made had two test cases, which are the Sandia Flame D (SFD) and the Verissimo Flameless Combustion Burner (FC). 
To better accommodate each test case some of the files were split into two versions:

1) EmissionsCalculator_Batch was split into three files: EmissionsCalculator_Batch (for the SFD), EmissionsCalculator_Batch_FC (for the FC) and EmissionsCalculator_Batch_plot (for plotting of the SFD).
   The main differences in the files are the test case boundary inputs and data dictionaries. The EmissionsCalculator_Batch_plot file is a more specific version of EmissionsCalculator_Batch.

2) PostProc was split into two files: PostProc (for the SFD) and PostProc_orig (for the FC).
   The main difference in the files is that in PostProc extra functions were added for plotting the SFD and in PostProc_orig one type of plot was added.

Furthermore there are a few things to note:

1) The Read_plot_all file was added that can read experimental and CFD data from the SFD.

2) The PostProcEmiss_cells was only used for the FC and can read the excel data for that case.
