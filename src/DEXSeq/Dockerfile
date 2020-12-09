FROM rocker/r-base:3.6.3 
LABEL description="Image for DEXSeq" 
LABEL maintainer="Alexander Dietrich" 


RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim python3 pip samtools
RUN R -e "install.packages(c('BiocManager','optparse'),dependencies=TRUE,repos='http://cran.rstudio.com/')"  
RUN R -e "BiocManager::install(c('DEXSeq','BiocParallel'),update=TRUE,ask=FALSE)" 

RUN pip install HTSeq pysam

WORKDIR /docker_main 


COPY [--chown=nobody:nogroup] ./DEXSeq.R /docker_main/DEXSeq.R 
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh 
RUN chmod 777 /docker_main/ENTRYPOINT.sh 
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]


