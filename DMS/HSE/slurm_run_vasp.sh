#!/bin/sh
##SBATCH --constraint=UBHPC&CPU-Gold-6130&INTEL&u25&OPA&MRI #request 'skylake' nodes
#SBATCH --constraint=MRI|NIH #request 'cascade' and 'skylake' nodes
#SBATCH --time=72:00:00
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
##SBATCH --qos=debug
##SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
##SBATCH --constraint=IB
#SBATCH --mem=128000
##SBATCH --mem=23000
#SBATCH --job-name=AgSSZ13-DMS-HSE
#SBATCH --output=job.out
#SBATCH --error=error.out
#SBATCH --mail-user=ruthbell@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

echo ">>>  VASP 5.4.4 with BEEF-vdw and vtst"
echo "--------------------------------------------------------"
echo ">>>  Job information"
echo "     SLURM_JOBID="$SLURM_JOBID
echo "     SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "     SLURM_NNODES"=$SLURM_NNODES
echo "     SLURMTMPDIR="$SLURMTMPDIR
echo "     working directory = "$SLURM_SUBMIT_DIR
echo "--------------------------------------------------------"
echo ">>>  Module information"
#Entering these module commands into the terminal prior to job submission removed errors:
#module load lmod
module load intel
module load mkl
#module load StdEnv
module load intel-mpi
############
module list
ulimit -s unlimited
#
echo "--------------------------------------------------------"
echo ">>>  Running Configuration"
# The initial srun will trigger the SLURM prologue on the compute nodes.
NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo NPROCS=$NPROCS
#The PMI library is necessary for srun
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
echo "--------------------------------------------------------"
echo ">>>  Job started at `date`"
echo "--------------------------------------------------------"
srun --propagate=STACK /projects/academic/mdupuis2/software/VASP_vtst_BEEF_vdw/vasp.5.4.4/build/std/vasp
echo "--------------------------------------------------------"
echo ">>>  Job finished at `date`"
