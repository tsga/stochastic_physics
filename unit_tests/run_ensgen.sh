#!/bin/sh
#SBATCH -e errlog
#SBATCH -o errlog
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=10
#SBATCH --job-name="ensgen"

module purge
source /home/Tseganeh.Gichamo/.my_mods

export OMP_NUM_THREADS=1
 
EnsForcGen=/scratch2/BMC/gsienkf/Tseganeh.Gichamo/stochastic_physics_mod/GenEnsForc.x

time srun '--export=ALL' --label -K -n 6 $EnsForcGen  
