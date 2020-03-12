#!/bin/bash

# Creates an index with every file whose name is *genome*.fa or *genome*.fasta

if [ -z "$abs_top_builddir" ]; then
    abs_top_builddir=..
fi

if [ -z "$abs_builddir" ]; then
    abs_builddir=.
fi

for file in $(find . -name "*genome*.fa" -o -name "*genome*.fasta"); do
    $LAUNCHER $abs_top_builddir/src/crac-index index $abs_builddir/${file%.*} $file
done
