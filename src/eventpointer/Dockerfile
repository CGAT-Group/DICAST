FROM rocker/r-base:3.6.3 
LABEL description="Image for EventPointer" 
LABEL maintainer="Alexander Dietrich" 


RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools


RUN R -e "install.packages(c('BiocManager','optparse'),dependencies=TRUE,repos='http://cran.rstudio.com/')"  
RUN R -e "BiocManager::install('EventPointer',update=TRUE,ask=FALSE)" 

RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y \
 python3.8 python3-pip

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./EventPointer.R /docker_main/EventPointer.R
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

