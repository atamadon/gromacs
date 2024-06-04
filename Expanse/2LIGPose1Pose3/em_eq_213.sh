#!/bin/bash
#SBATCH -J em_eq_2lig_pose1pose3                 # Job name
#SBATCH -o em_eq.o%j                # Name of stdout output file
#SBATCH -e em_eq.e%j                # Name of stderr error file
#SBATCH -p compute                  # Queue (partition) name
#SBATCH -A ucb344                   # Account
#SBATCH --export=ALL
#SBATCH --nodes=5                   # Number of nodes   
#SBATCH --ntasks-per-node=8         # Number of tasks per node
#SBATCH --cpus-per-task=8
#SBATCH --mem=249000M
#SBATCH -t 03:00:00                 # This job will run for a maximum of 24 hrs
#SBATCH --mail-user=atamadon@berkeley.edu
#SBATCH --mail-type=all             # Send email at begin and end of job

# This job runs with 5 nodes, 8 cores per node for a total of 40 cores.
# We use 8 MPI tasks and 8 OpenMP threads per MPI task (gmx may limit to 64 tasks)

# # Loading necessary modules
# module purge
# module load slurm
# module load cpu/0.17.3b
# module load aocc/3.2.0/io3s466
# module load openmpi/4.1.3/xigazqd

# # Loading Gromacs version
# module load gromacs/2020.4/2ufeq67-omp

# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the GROMACS simulation

EM_MDP="../Sample_mdp_files/minim.mdp"
NVT_MDP="../Sample_mdp_files/nvt1.mdp"
NPT_MDP="../Sample_mdp_files/npt1.mdp"
MD_MDP="../Sample_mdp_files/production.mdp"

# Energy Minimization
# Local
gmx grompp -f $EM_MDP -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
# HPC-----------------------------------------
# gmx_mpi grompp -f $EM_MDP -c solv_ions.gro -r solv_ions.gro -p topol.top -o em.tpr -maxwarn 2
# mpirun -np 32 gmx_mpi mdrun -v -deffnm em
# echo 13 0 | gmx_mpi energy -f em.edr -o em_Epot.xvg

# NVT Equilibration
# Local---------------------------------------
gmx grompp -f $NVT_MDP -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -v -deffnm nvt
# HPC-----------------------------------------
# gmx_mpi grompp -f $NVT_MDP -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2
# mpirun -np 40 gmx_mpi mdrun -v -deffnm nvt

# NPT Equilibration
# Local---------------------------------------
gmx grompp -f $NPT_MDP -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -v -deffnm npt
# HPC-----------------------------------------
# gmx_mpi grompp -f $NPT_MDP -c nvt.gro -r nvt.gro -t nvt.cpt -n index.ndx -p topol.top -o npt.tpr -maxwarn 2
# mpirun -np 40 gmx_mpi mdrun -v -deffnm npt
