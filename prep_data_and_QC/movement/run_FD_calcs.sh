#!/usr/bin/env bash

# Load matlab
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Matlab2019b --silent

matlab2019b -nodisplay -nodesktop -r "run UCLA_CNP_calc_FD.m, exit"
# matlab2019b -nodisplay -nodesktop -r "run ABIDE_ASD_calc_FD.m, exit"