#!/bin/sh
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t 2:00:00
#SBATCH -J buta_tddft
#SBATCH -o buta_tddft.txt
#

export PATH=$PATH:/software/orca/4.0.1.2
#export RSH_COMMAND="/usr/bin/ssh -x"


cp buta_tddft.inp /scratch/zpiontko/buta/tddft/
cd /scratch/zpiontko/buta/tddft/

/software/orca/4.0.1.2/orca buta_tddft.inp > buta_tddft.out

cp * /home/zpiontko/data/Spring_2017/20171129_deltex_project/butadiene_orcaRRcalcs
