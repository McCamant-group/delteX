#!/bin/sh
#SBATCH -p standard
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t 2:00:00
#SBATCH -J buta_opt_freq
#SBATCH -o buta_opt_freq.txt
#

export PATH=$PATH:/software/orca/4.0.1.2
#export RSH_COMMAND="/usr/bin/ssh -x"


cp buta_opt_freq.inp /scratch/zpiontko/buta/opt_freq/
cd /scratch/zpiontko/buta/opt_freq/

/software/orca/4.0.1.2/orca buta_opt_freq.inp > buta_opt_freq.out

cp * /home/zpiontko/data/Spring_2017/20171129_deltex_project/butadiene_orcaRRcalcs
