# Booster_Shot_Variant
This repository contains original code created for the study titled "Antigen presentation dynamics shapes the response to emergent variants like the Omicron strain of SARS-CoV-2 after multiple vaccinations with the wild type strain" by Yang, Van Beek, Wang, Muecksch, Canis, Hatziioannou, Bieniasz, Nussenzweig, and Chakraborty in Cell Reports (2023). [Link to the paper](https://www.cell.com/cell-reports/fulltext/S2211-1247(23)00267-X) The code has also been published in Mendeley Data ([Link](https://data.mendeley.com/datasets/39bb2273yz))

The raw data from simulations can be downloaded from the following [link](https://mitprod-my.sharepoint.com/:u:/g/personal/leerang_mit_edu/EdtAr52PD-5IoVeTmPcYGXEBVk8vmUmOnXmwJ6IIxzIL_A?e=2DlfJh ). The dataset is large (111 GB).

The results from this study can be reproduced through the following steps.  
1. Run “Code_Parameter_Generation/getParameters.m” to generate text files that contain the parameters that are given as inputs to the simulations. These text files will be saved in the folder “parameters”. Each text file is used to run a group of simulations that are related. Each row of a text file gives a set of parameters for a single simulation.  
2. Use these parameters as inputs to run “Code_Simulation/runGCs.m”, or run “Code_Simulation/runGCs_All.m”. The former should be used for text files that start with a specific dose number (e.g. “Vax1..”). The latter should be used for other text files. To run simulations with all sets of parameters, using a parallel computing cluster will be necessary. If the cluster uses the Slurm job scheduling system, then the following detailed instruction can be followed. Otherwise, appropriate steps should be taken to run the simulations.  

<blockquote>
If the cluster uses Slurm:  <br/>
  
(1)	Run “Code_Parameter_Generation/getJobsToRun.m”. This will generate a text file “jobs_to_run.txt” in the current folder. This file contains lines of Linux commands that can be executed to submit batches of jobs to Slurm.  
(2)	Using the above file, submit all the jobs that correspond to Vax1. When all these jobs finish, which will take about 1-2 hours if all jobs run parallelly, then submit the jobs that correspond to Vax2. These jobs will take about 5-10 hours to finish if all jobs run parallelly. Afterward, submit the jobs that correspond to Vax3. Finally submit the jobs that correspond to Vax4.  
Once the simulations are all run, the results from the simulations are saved in the folder “Data”, under appropriate subfolders that describe the simulation conditions, in the form of MATLAB data files (.mat files).  
</blockquote>

3. The codes that analyze and plot results are in the folder “Code_Plotting”. It contains subfolders that correspond to each main and supplemental figure. In these subfolders, there are two scripts for each sub panel. The script with a longer name (e.g. “Figure2b_GCOutcomes.m”) analyzes the data and then saves a summary “.mat” file in “Code_Plotting/Results” folder with the same name as the script file. This file contains all the information needed to plot the specific subpanel. Then, the script with a shorter name (e.g. “Figure2b”) reads this file and plots the panel.

* Equation 8 of STAR Methods shows that the effects of mutations on the left-hand side are proportional to random variables sampled from shifted log-normal distribution. To obtain equality relationship, the right-hand side should be multiplied by -log10(e). See the code in Code_Simulation/getNaive for the implementation.
