#!/bin/bash
ImageArray=( "$@" )
slurmarraylength=${#array[@]}

cat<<EOF>./slurmsubmit.sh
#!/bin/bash
#SBATCH --array=1-${slurmarraylength}%1
EOF

# Get SBATCH directives from slurm-docker-runner.sh
    let cnt=$(cat ./slurm-docker-runner.sh |wc -l)-1
cat slurm-docker-runner.sh |tail -${cnt} |head -n6 >> ./slurmsubmit.sh

# Make a list of images to run and select specific image to run
echo 'image=$(echo '${ImageArray[@]}' | cut -d" "  -f${SLURM_ARRAY_TASK_ID})' >> ./slurmsubmit.sh

#Final echo & docker run commands to complete slurmsubmit.sh
tail ./slurm-docker-runner.sh -n2 >> ./slurmsubmit.sh