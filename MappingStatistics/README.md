# Mapping Statistics
The scripts take a gtf file, to parse the transcriptomic read coordinates of an AS-Simulator output to genomic coordinates.
Then these reference reads are compared to the mapping in a .bam file.

## Script usage with console
Rscript main_server.R

All 4 files (main_server.R, config.txt, functions.R, read_fastq.awk) have to be in the same directory.

Input files and output locations get specified in config.txt.

The variables in config.txt are:

bam_path: A path to a bam file, containing the mapped reads.

fastq_path: A file to a fastq file with the reference reads. The headers have to be in AS-Simulator output format.

gtf_path: The gtf file used by the AS simulator for read simulation.

genomic_reads: A .RData file containing the reference reads with genomic coordinates. If it doesn't already exist, it will be created.
If the file does exist, its data will be loaded, to skip the coordinate parsing process.

output_folder: The output directory name. (Directory has to be in the same folder as the scripts. If it doesn't exist, it will be created.)

## Output

The output folder will contain 3 files:

bam_reads.RData: The bam file reads in genomic ranges format

full_info.csv: Mapping quality information for each individual read

summary.csv: Summary of full_info.csv


summary.csv contains the following measurements:

n_reads: The number of reference reads

n_proper_maps: The number of mapped reads

n_correct_chrom: The number of reads, mapped to the right chromosome

n_correct_start: The number of reads, with correct start position

n_correct_end: The number or reads, with correct end position

n_correct_range: The number or reads, with correct start and end position

n_all_junctions_correct: The number of correctly mapped reads

n_junctions: The number of all junctions of the reference reads

n_correct_junctions: The number of correctly mapped junctions

n_false_junctions: The number of falsely mapped junctions

n_pairs: The number of reference read pairs

n_proper_pairs: The number of mapped read pairs

n_correct_chrom_pairs: The number of read pairs, mapped to the correct chromosome

n_correct_range_pairs: The number of read pairs, where both read ranges are correct

n_all_junctions_correct_pairs: The number of correctly mapped read pairs

## Plotting
The R-notebook plots_mapping.Rmd can be used to visualize the results.

## Required libraries
The scripts utilize the following R libraries:
- Rsamtools
- data.table
- GenomicFeatures
- tidyr
- GenomicAlignments
- pryr
- ggplot2


