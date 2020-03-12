#!/bin/bash

# Test that the indexes genomes are identical to the ones retrieved using the
# get option from crac-index

if [ -z "$abs_top_builddir" ]; then
    abs_top_builddir=..
fi

if [ -z "$abs_builddir" ]; then
    abs_builddir=.
fi

compare_dna_seq() {
    raw_seq1=$(tempfile)
    raw_seq2=$(tempfile)
    (grep -v '^>' $1 | tr -dc 'ACGTTacgtn ' | tr 'acgtn ' 'ACGTNN'; echo) > $raw_seq1
    (grep -v '^>' $2 | tr -dc 'ACGTTacgtn ' | tr 'acgtn ' 'ACGTNN'; echo) > $raw_seq2
    diff=$(cmp -b -l $raw_seq1 $raw_seq2 | awk 'BEGIN{c=0} $3 != "N" && $5 != "N"{c++} END {print c}')
    rm -f $raw_seq1 $raw_seq2
    echo $diff
}

error=0
TEMP=$(tempfile)
nb_files=$(find $abs_builddir -name "*.ssa" | wc -l)
i=1
error=0
{
echo "1..$nb_files"
for file in $(find $abs_builddir -name "*.ssa"); do
    file=${file%.ssa}
    $LAUNCHER $abs_top_builddir/src/crac-index get $TEMP $file
    orig_fasta="."${file#$abs_builddir}
    orig_fasta=${orig_fasta}.fa*
    if [ -f $orig_fasta ]; then
        orig_fasta=$(ls $orig_fasta)
        if [ $(compare_dna_seq $TEMP $orig_fasta) -gt 0  ]; then
            echo -n "not "
        fi
    else
        echo -n "not "
        error=$((error+1))
    fi
    echo "ok $i - file "$orig_fasta
    i=$((i+1))
done
} > ${0%.sh}".tap"

rm -f $TEMP

[ $error -eq 0 ]
