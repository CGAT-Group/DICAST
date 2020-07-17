FROM rocker/r-base:4.0.0
LABEL description="Image for SGSeq"
LABEL maintainer="Alexander Dietrich"

RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim

RUN R -e "install.packages(c('BiocManager','optparse'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('SGSeq',update=TRUE,ask=FALSE)" 

WORKDIR /docker_main

ENTRYPOINT ["/bin/bash"]
