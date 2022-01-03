#!/bin/bash

# This is an example of the inputs required to run DICAST on the Human Genome.

# Building the directory structure needed for DICAST.
source scripts/config.sh
mkdir -p $(echo $controlbam| sed 's/^\/MOUNT/./g') $(echo $controlfastq| sed 's/^\/MOUNT/./g') $(echo $bowtie_fastadir| sed 's/^\/MOUNT/./g')


# Downloading Human references fron Ensemble's ftp.
link="http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/"
## Downloading bowtie genome fastas for each Chromosome.
for chromosomes in $(curl $link | cut -d ' ' -f2 | cut -d '"' -f3 | grep -v "nonchromosomal\|primary\|toplevel\|dna_\|alt" | grep Homo_sapiens|sed 's/...>//g'| tr -d '>'); do curl -o $(echo $bowtie_fastadir| sed 's/^\/MOUNT/./g')/$(echo $chromosomes | cut -d '.' -f 5-) $link$chromosomes; done

## Downloading full genome as 1 fasta.
curl -o input/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#Downloading gff files.
curl -o input/Homo_sapiens.GRCh38.105.gtf.gz http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz
curl -o input/Homo_sapiens.GRCh38.105.gff3.gz http://ftp.ensembl.org/pub/release-105/gff3/homo_sapiens/Homo_sapiens.GRCh38.105.gff3.gz