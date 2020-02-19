<img src="https://d1q6f0aelx0por.cloudfront.net/product-logos/library-docker-logo.png" alt="Docker">

This is a list of Dockers made for the Benchmarking pipeline that [exbio](https://www.baumbachlab.net/) will build for Sys_Ca_Re. These dockers would also be a valuable resource for all to use, in Bioinformatics. 


##Usage
These images should all be usable by the command:
```shell
cd dockers/*/

docker build ./ -t amit/**TOOLNAME**

cd path/to/fastq/files/

docker run -v "$(pwd)":/myvol1 --user $(id -u):$(id -g)  amit/

```
