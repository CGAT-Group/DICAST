#!/bin/bash

tool=psisigma
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

mk_outdir

create_input(){
	folder=$1
	output_file=$2
  if ! [[ -f output_file ]]; then
    find $folder/ -maxdepth 1 -type f -name '*.bam' >> "$output_file"
  fi
}

cd afolder
ln -s $renamed_casebam/*.bam* .
ln -s $renamed_controlbam/*.bam* .
ln -s $renamed_star_alignment_files/*.SJ.* .

create_input $renamed_casebam groupa.txt
create_input $renamed_controlbam groupb.txt

sorted_gtfname="sorted_${gtfname}"

if [[ ! -f "$sorted_gtfname" ]]; then
  echo "Sorting gtf file $gtfname..."
  (grep "^#" "$gtf"; grep -v "^#" "$gtf" | sort -k1,1 -k4,4n) > "$sorted_gtfname"
  echo "Done sorting."
fi

psi_sigma_script=../PSI_Sigma/dummyai.pl
perl $psi_sigma_script --gtf $sorted_gtfname --name PSIsigma --type 1 -nread 10

cp PSIsigma_r10_ir3* $outdir
cp *.IR.out.tab $outdir
cp *.db $outdir
