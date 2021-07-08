#!/bin/bash
tools=(asgal aspli eventpointer irfinder majiq sgseq spladder whippet )
pushd $1
for j in ${tools[@]} 
do 
	echo $j && \
		for i in $(find ./ -name '*unified.out')
		do basename $i
		done | grep $j |wc -l \
			&& for i in $(find ./ -name '*unified.out')
			do basename $i
			done | grep $j
	echo ------------
done

for j in ${tools[@]} 
do 
	echo $j && \
		for i in $(find ./ -name '*comparison.txt')
		do basename $i
		done | grep $j |wc -l \
			&& for i in $(find ./ -name '*comparison.txt')
			do basename $i
			done | grep $j
	echo -----------
done
popd