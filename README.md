# CMIP6_data
 Efficient access for CMIP6 data structures.

Here is a bunch of python code for accessing CMIP6 data on the jasmin research cluster.
It may also work well on other clusters, you'd have to check how the CMIP6 data directory structure works.

The files are as follows: 
CMIP_functions.py  
CMIP_inputs.py  
grid_set.py  
Essential modules for finding files and reading them. There's lots of regridding too.

Process_CMIP6.py
The example processing file. Lots of annotation in here

process_single_job
sbatch script for using the jasmin queuing system, runs a single job
use 'sbatch process_single_job'

batch_process_jobs
sbatch script for using the jasmin queuing system, runs a whole batch of jobs (good for ensembles)
use 'sbatch batch_process_jobs'

proc_inputs/proc_input1.txt
proc_inputs/proc_input2.txt
two example input files for Process_CMIP6.py. They contain information for two working scenarios

run_outs/
batch_run_outs/
these contain all the Processing output from the two example runs

Outputs/
two example NetCDF files containing the example runs

CMIP6_open_processed.ipynb
A notebook that will open and visualise the processed files.


add_scenarios.py  
This is an extra, undocumented data downloading script. it will build a mirror of the jasmin CMIP6 archive, and allow you to download extra data to that found on Jasmin. The existing data is sym linked, extra data is added. Then new structure can be accessed by the processing script as before.
