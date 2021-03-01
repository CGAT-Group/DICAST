<img src="https://gitlab.lrz.de/ge46ban/dockers/-/blob/develop/docs/source/img/icon.png" alt="DICAST">

This is a list of Dockers made for the Benchmarking pipeline that [exbio](https://www.baumbachlab.net/) will build for Sys_Ca_Re. These dockers would also be a valuable resource for all to use, in Bioinformatics. 

### Prerequisites

These dockers work on input files in a working directory. 

#### Mapping tools

Input files: fastq files 



Other required files:  fasta files, unzipped: "*.fa" and annotation files, unzipped: "*.gtf"

### Usage
These images should all be usable by the command:
```shell
cd dockers/*/

docker build ./ -t amit/**TOOLNAME**

cd path/to/fastq/files/

docker run -v "$(pwd)":/myvol1 --user $(id -u):$(id -g)  amit/

```

If you'd like a cheat sheet on using dockers with this repo, it's available [here](https://gitlab.lrz.de/ge46ban/dockers/-/wikis/Docker-commands-CHEAT-SHEET)
