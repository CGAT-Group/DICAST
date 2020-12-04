FROM rocker/r-base:3.6.3
LABEL description="Image for leafcutter"
LABEL maintainer="Alexander Dietrich"

WORKDIR /opt

RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y git samtools curl libcurl4-openssl-dev libssl-dev libxml2-dev vim software-properties-common g++ gsl-bin libgsl-dev libv8-dev cmake
RUN apt-get --no-install-recommends --fix-broken install -y python2.7 python2.7-dev
RUN apt-get --no-install-recommends --fix-broken install -y r-cran-rstan

RUN git clone https://github.com/griffithlab/regtools && cd regtools/ && mkdir build && cd build/ && cmake .. && make

WORKDIR  /docker_main/

#install leafcutter
RUN R -e "install.packages(c('devtools','BiocManager','V8'))"
RUN R -e "BiocManager::install(c('Biobase','DirichletMultinomial'))"
RUN R -e "devtools::install_github('stan-dev/rstantools')"

RUN R -e "devtools::install_github('davidaknowles/leafcutter/leafcutter')"

#download git-repo for python scripts
RUN git clone https://github.com/davidaknowles/leafcutter


ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/regtools/build

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

