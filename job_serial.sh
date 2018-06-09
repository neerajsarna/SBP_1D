#!/usr/bin/env zsh
 
### Job name
#BSUB -J "MATLAB_ARRAY[26,28,30,32,34,36,38,40,27,29,31,33,35,37,39,41]"
 
### File / path where STDOUT will be written, the %J is the job id
#BSUB -o log_files/solving_gaussian_collision_%I
 
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
matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile log_files/solve_collision_gaussian_1x1v_$LSB_JOBINDEX.log <<EOF
run solve_collision_gaussian_1x1v($LSB_JOBINDEX);
quit();
EOF
