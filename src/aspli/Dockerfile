FROM rocker/r-base:3.6.3
LABEL description="Image for ASpli"
LABEL maintainer="Alexander Dietrich"


RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim

RUN R -e "install.packages(c('optparse','BiocManager'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('rtracklayer',update=TRUE,ask=FALSE)"
RUN R -e "BiocManager::install('ASpli',update=TRUE,ask=FALSE)"

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./ASpli.R /docker_main/ASpli.R
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]


