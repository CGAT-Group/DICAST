FROM rocker/r-base:3.6.3
LABEL description="Image for KisSplice and KissDE"
LABEL maintainer="Alexander Dietrich"


#############
# KisSplice #
#############

#install python3
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y python3 python3-dev python3-pip

#download and install KisSplice in /opt
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y cmake zlib1g
WORKDIR /opt
RUN wget ftp://pbil.univ-lyon1.fr/pub/logiciel/kissplice/download/kissplice-2.5.1.tar.gz && tar -xzf /opt/kissplice-2.5.1.tar.gz && cd /opt/kissplice-2.5.1 && cmake . && make

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/kissplice-2.5.1/bin


#############
#  KissDE   #
#############

RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev 

RUN R -e "install.packages('BiocManager',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('kissDE',update=TRUE,ask=FALSE)" 

WORKDIR /docker_main

ENTRYPOINT ["/bin/bash"]

