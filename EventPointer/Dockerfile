FROM rocker/r-base:3.6.3 
LABEL description="Image for ASpli" 
LABEL maintainer="Alexander Dietrich" 


RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev 


RUN R -e "install.packages('BiocManager',dependencies=TRUE,repos='http://cran.rstudio.com/')" 
#RUN R -e "BiocManager::install('rtracklayer',update=TRUE,ask=FALSE)" 
RUN R -e "BiocManager::install('EventPointer',update=TRUE,ask=FALSE)" 

ENTRYPOINT ["/bin/bash"]
