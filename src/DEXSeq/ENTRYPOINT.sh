#!/bin/bash 

tool=dexseq 
source /MOUNT/scripts/config.sh 
source /MOUNT/scripts/asevent_config.sh 
source /MOUNT/scripts/asevent_func.sh



test_gtf $gtf

#make output directory
mk_outdir $tool

#handle SAM files
handlesamfiles 1

#use dexseq python script to create gff version of gtf file with same name as gtf
python_script_paths=/usr/local/lib/R/site-library/DEXSeq/python_scripts
gffname=$(basename $gtf .gtf).gff

echo preparing annotation...
python3 $python_script_paths/dexseq_prepare_annotation.py $gtf $outdir/$gffname


#count reads for bamfiles in case-folder
mkdir $outdir/tmp
chmod -R 777 $outdir/tmp

echo counting reads in case-folder...
readbamfiles $casebam casebamlist
for filename in $(cat $outdir/casebamlist)
do
	out_name=$(basename $filename .bam).g1.txt
	python3 $python_script_paths/dexseq_count.py -f bam $outdir/$gffname $filename $outdir/tmp/$out_name
done

#count read for bamfiles in control-folder
echo counting read in control-folder
readbamfiles $controlbam controlbamlist
for filename in $(cat $outdir/controlbamlist)
do
        out_name=$(basename $filename .bam).g2.txt
        python3 $python_script_paths/dexseq_count.py -f bam $outdir/$gffname $filename $outdir/tmp/$out_name
done


echo starting DEXSeq differential exon usage analysis...
Rscript /docker_main/DEXSeq.R --gff $outdir/$gffname --out $outdir --ncores $ncores

cleaner

