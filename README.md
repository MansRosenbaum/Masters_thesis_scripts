# Masters_thesis_scripts

This repository includes the scripts I have made both to automize and run MolDStruct, but also for creating structures from E-feild simulation data, doing analyses, fail-searching and certain plots. 

The data from simulations of proteins exposed to an external electric field was taken from A.Sinelnikova, https://doi.org/10.1039/D0SC06008A. This includes E-field simulation trajectory files and topology files.

Separate codes were created to automize the MolDStruct pipeline depending on what datasets that were used. The same was done for analysis.

**Descriptions**

analysis_all_datasets.ipynb
 - Notebook for creating ion maps, movies and using dimensionality reduction methods for all-dataset runs.                                           

analysis_correlations_all_datasets.ipynb
 - Notebook for creating ion maps, movies and finding correlations described in the report.                       

analysis_single_dataset.ipynb
 - Notebook for creating ion maps, movies and using dimensionality reduction methods for single-dataset runs.

analysis_2runs.ipynb
- An analysis-notebook which shows the effect of putting two single datasets (Run 0 and Run 8) in the same space before using dimensionality reduction techniques. These results are not used in the report.        

failed_runs.ipynb
 - Script made to search for the MolDStruct folders for failed simulations. This could tell the user how many fails the were, what folder they are in and were the fail occured (by comparing the content of the folder with a successful run, revealing what files that are missing in the failed run). These were only a small fraction of the total number of simulations and are thus not significant in this thesis.

gaussian_pulse.ipynb
 - Notebook for plotting a Gaussian curve symbolizing the laser pulse in MolDStruct, Figure 5 in the report.

make_structures.ipynb
 - Notebook for exctracting the 101 structures from E-field simulation trajectory files.

mds_run_all_datasets.py
 - Python script that runs MolDStruct on all 101 structure files of all 100 E-field simulation runs, making a total of 10100 MolDStruct simulations.

mds_run_single_dataset.py
 - Python script that runs MolDStruct on all 101 structure files for one E-field run with an adjustable number of simulations per timestep.

plot_cross_sections.ipynb
 - Notebook for plotting cross sections for photo-absorption at a 0-1500 eV interval. Cross-sections are taken from https://vuo.elettra.eu/services/elements/WebElements.html.

run_moldstruct.py
 - Python script for running single MolDStruct simulations. The code can automatically download pdb-files given a PDB-code (such as 1UBQ). Local pdb- or gro-files can also be implemented.

**Parameters**

Parameters for running MolDStruct simulations such as timesteps, photon energy, the number of photons and FWHM can be adjusted in the scripts. Helper functions will then enter the mdp-file and change this automatically. run_moldsctruct.py also has interactive commands for the user to implement preferable parameters when running the code. For convenience when running the code multiple times, this was removed in mds_run_single_dataset.py and mds_run_all_datasets.py. In both of these scripts, you can however parallellize multiple runs. In these codes you can change the number of CPU cores you want to utiilize, as well as the number of simulations for each timestep in mds_run_single_dataset.py. In mds_run_all_datasets.py, the number of simulations per structure is set to 1, making 10100 simulations in total.

You must have MolDStruct downloaded beforehand and change the directory pathway.





