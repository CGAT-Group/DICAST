FROM  dicast/star:0.2 as builder

FROM rocker/r-base:4.0.0
LABEL description="Image for SGSeq"
LABEL maintainer="Alexander Dietrich"

ARG star_version=2.7.5c
#Copying STAR from builder docker
COPY --from=builder /docker_main/STAR-${star_version} /opt/STAR-${star_version}/
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/STAR-${star_version}/bin/Linux_x86_64

RUN R -e "install.packages(c('BiocManager','optparse'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('SGSeq',update=TRUE,ask=FALSE)" 

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./SGSeq.R /docker_main/SGSeq.R
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]