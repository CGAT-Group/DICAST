FROM rocker/r-base:4.0.0
LABEL description="Image for SGSeq"
LABEL maintainer="Alexander Dietrich"

ARG star_version=2.7.3a
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean

RUN R -e "install.packages(c('BiocManager','optparse'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('SGSeq',update=TRUE,ask=FALSE)" 

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

