#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=dockermappers
#SBATCH --output=/nfs/home/users/afenn/slurm/logs/%x.%j.%a.out
#SBATCH --error=/nfs/home/users/afenn/slurm/logs/%x.%j.%a.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=5

# This script needs an array $images with a list of images made by arraymaker.sh
# This scripts needs to assign $image done with arraymaker.sh
# Do not modify this code, line numbers are counted by other scripts.

echo $image $USER
docker run -v "/nfs/proj/AS_dockers/:/MOUNT/" --user $(id -u):$(id -g) --rm proj/${image}:0.01 ; echo This docker exited with status is: $?