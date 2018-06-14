#!/usr/bin/env zsh
 
### Job name
#BSUB -J "MATLAB_ARRAY[43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,44,46,48,50,52,54,56,58,60,62,64,66,68,70]"
 
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
