#!/bin/bash
#SBATCH --array=1-22%1
#SBATCH --ntasks=1
#SBATCH --job-name=dockermappers
#SBATCH --output=/nfs/home/students/fanny/slurm/logs/%x.%j.%a.out
#SBATCH --error=/nfs/home/students/fanny/slurm/logs/%x.%j.%a.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=5

# Making an array $images with a list of images

image=${1:-$(cat /nfs/proj/AS_dockers/images.txt | head -${SLURM_ARRAY_TASK_ID}|tail -1 )}
echo $image $USER
docker run -v /nfs/proj/AS_dockers/unify_test_MOUNT/:/MOUNT --user $(id -u):$(id -g) --rm $image ; echo This docker exited with status: $?
