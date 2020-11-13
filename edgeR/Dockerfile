FROM rocker/r-base:4.0.0
LABEL description="Image for edgeR"
LABEL maintainer="Alexander Dietrich"

RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools

RUN R -e "install.packages(c('optparse','BiocManager','statmod'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('edgeR','Rsubread'),update=TRUE,ask=FALSE)"

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./edgeR.R /docker_main/edgeR.R
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

