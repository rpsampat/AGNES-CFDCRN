This code was written by Rishikesh Sampat as part of a MSc thesis in Aerospace Engineering at TU Delft in association with Mitsubishi Turbocharger and Engine Europe B.V.

The code is written in Python 2.7. The following Python files are part of the package:
1) casefilepy
2) datfilepy
3) CRN_Gen
4) BFS_pchanged
5) ReactorGen_min_internalmfc_outletvalve_newton
6) SaveReactors
7) PostProc
8) PostProcEmiss_cells
9) EmissionsCalculator_Batch

1) casefilepy: Reads the ANSYS Fluent case file that needs to be in non-binary format. It is possible to convert a binary case file to non binary by simply unselecting the 'binary' checkbox in the 'SaveAs' option in ANSYS Fluent.

2) datfilepy: Reads the ANSYS Fluent data file that needs to be in no-binary format.

3) CRN_Gen: Reads in the case and data graphs of the mesh and prepares for clustering. The clustering code, BFS-Pchanged, is called iteratively for the clustering process. The final clustered graph is stored with the corresponding data.

4) BFS_pchanged: Clustering algorithm implemented.

5) ReactorGen_min_internalmfc_outletvalve_newton: CRN solver.

6) SaveReactors: Saves reactor data in Excel file and creates 3d CRN representation in ParaView.

7) PostProc: Creates a 2D plane of CRN data projected back to the CFD domain.

8) PostProcEmiss_cells: Python plots at different locations in the combustion chamber.

9) EmissionsCalculator_Batch: Master file to run the CFD-CRN process.

-----------------------------------------------------------------------------
The ANSYS FLuent case and data files should be present in the current folder with the name 'CFD_CRN'. The order of executing the program is :
1)casefilepy
2)datfilepy
3)CRN_Gen
4)ReactorGen_min_internalmfc_outletvalve_newton

This is the same order as executed by EmissionsCalculator_Batch.py and the user may directly use this program.
When different CRNs need to be generated from the same CFD simulation, 1) and 2) need to be done only for the first case.

Executing the code creates a subfolder within the current folder containing the data of the created CRN. The folder has the following naming convention:
'number of reactors specified_number of reactors created_clustering criteria'.

For example:
'6000_6882_Static Temperature_Mass fraction of h2o_Mass fraction of ch4_Velocity Dir_doe_reactnum'.

Within the subfolder several files are created. The final output data files are named 'Emission10.xls' and 'data10.xls'.
1)data10: data of all the reactors, species, temperature etc.
2)Emissions10: NOx at certain locations in the system. The values are uncorrected.
3)NumberOfReactors1: Number of reactors at the end of each iteration of clustering
4)ReactorCells1: Number of cells clustered in each reactor
5)Tolerance1: Tolerance of clustering at the end of each clustering iteration.