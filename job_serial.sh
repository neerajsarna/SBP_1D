#!/usr/bin/env zsh
 
### Job name
#BSUB -J "MATLAB_ARRAY[3,7,11,15,19,23,27,4,8,12,16,20,24,44,45]"
 
### File / path where STDOUT will be written, the %J is the job id
#BSUB -o log_files/solving_inflow_%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 5:00
 
### Request memory you need for your job in MB
#BSUB -M 5000
 
### Change to the work directory
cd /home/ns179556/SBP_1D/
 
### load modules and execute
module load MISC
module load matlab
 
 
# start non-interactive batch job
matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile log_files/solve_inflow_1x3v_$LSB_JOBINDEX.log <<EOF
run solve_inflow_1x3v($LSB_JOBINDEX);
quit();
EOF
