#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --exclusive
#
#This recomend to minimize clashes with other users
#SBATCH --partition=scaling
#

cmd_line="sed -n '$SLURM_ARRAY_TASK_ID p' $1"
echo "get line: $cmd_line"

cmd=$(eval $cmd_line)
echo "command: $cmd"

run_cmd=$(eval "echo '$cmd' | cut -d \> -f 1")
stdout=$(eval "echo '$cmd' | cut -d \> -f 2")
echo "run command: $run_cmd"
echo "output: $stdout"

eval $cmd

#srun --output=$stdout $run_cmd


