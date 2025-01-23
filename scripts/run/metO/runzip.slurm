#!/bin/bash
#SBATCH --job-name=NAARC-ZIP
#SBATCH --time=24:00:00
#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

. XXX_ENV_XXX

if [ ! -d OUTPUTS_ZIP ]
then
	mkdir OUTPUTS_ZIP
fi

if [ ! -d OUTPUTS_PROCESSED ]
then
	mkdir OUTPUTS_PROCESSED
fi

cd OUTPUTS

for i in NAARC_1d_** NAARC_1m_*
do
        if [[ $fname == *"grid_T"* ]]; then
	   nccopy -4 -d 1 -c x/128,y/128,deptht/75 $i ../OUTPUTS_ZIP/$i
	elif [[ $fname == *"grid_U"* ]]; then
           nccopy -4 -d 1 -c x/128,y/128,depthu/75 $i ../OUTPUTS_ZIP/$i
	elif [[ $fname == *"grid_V"* ]]; then
           nccopy -4 -d 1 -c x/128,y/128,depthv/75 $i ../OUTPUTS_ZIP/$i
	elif [[ $fname == *"grid_W"* ]]; then
           nccopy -4 -d 1 -c x/128,y/128,depthw/75 $i ../OUTPUTS_ZIP/$i
	fi

	mv $i ../OUTPUTS_PROCESSED/

done

