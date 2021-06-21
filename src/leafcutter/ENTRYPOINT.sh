#!/bin/bash

# use config and function file
tool=leafcutter
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

lc_folder=/docker_main/leafcutter

#test input-files
test_gtf $gtf

#making outdir
mk_outdir 

#handle samfiles
handlesamfiles 1

trap cleaner EXIT

#convert bam to junc files:
#also create grouping txt file 
mkdir -p $outdir/junc_files
readbamfiles $casebam casebamlist
for filename in $(cat /tmp/casebamlist)
do 
	name=$(basename -s .bam $filename)
	touch $outdir/junc_files/$name.junc
	outfile=$outdir/junc_files/$name.junc
	regtools junctions extract -a 8 -s 0 -m 50 -M 500000 $filename -o $outfile
	echo $outfile >> $outdir/casejunc.txt
	echo -e "$name\tCASE" >> $outdir/groups_file.txt
done

readbamfiles $controlbam controlbamlist
for filename in $(cat /tmp/controlbamlist) 
do 
        name=$(basename -s .bam $filename)
	touch $outdir/junc_files/$name.junc
        outfile=$outdir/junc_files/$name.junc
	regtools junctions extract -a 8 -s 0 -m 50 -M 500000 $filename -o $outfile
        echo $outfile >> $outdir/controljunc.txt
	echo -e "$name\tCONTROL" >> $outdir/groups_file.txt
done

#combine junc.txt files
cat $outdir/*junc.txt > $outdir/junc_all.txt

mkdir -p $outdir/clustering
python2.7 $lc_folder/clustering/leafcutter_cluster_regtools.py -j $outdir/junc_all.txt -m 50 -r $outdir/clustering -o CASEvsCONTROL -l 500000 -k

#convert gtf to exon-file
#requires gtf in gzip format
gzip -c $gtf > $outdir/gtf.gz
Rscript $lc_folder/scripts/gtf_to_exons.R $outdir/gtf.gz $outdir/exon_file.txt


#differential splicing analysis
Rscript $lc_folder/scripts/leafcutter_ds.R  --num_threads $ncores --exon_file $outdir/exon_file.txt $outdir/clustering/CASEvsCONTROL_perind_numers.counts.gz $outdir/groups_file.txt -o $outdir


rm -f $outdir/exon_file.txt
rm -f $outdir/junc_all.txt
rm -f $outdir/groups_file.txt
rm -f $outdir/casejunc.txt
rm -f $outdir/controljunc.txt
rm -f /tmp/casebamlist
rm -f /tmp/controlbamlist
rm -f $outdir/gtf.gz