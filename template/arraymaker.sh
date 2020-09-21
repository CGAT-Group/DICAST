#!/bin/bash
arno=$( cat /nfs/proj/AS_dockers/images.txt |wc -l)
cat<<EOF>./slurmsubmit.sh
#!/bin/bash
#SBATCH --array=1-${arno}%1
EOF
let cnt=$(cat ./slurm-docker-runner.sh |wc -l)-1
cat slurm-docker-runner.sh |tail -$(echo $cnt ) >> ./slurmsubmit.sh

