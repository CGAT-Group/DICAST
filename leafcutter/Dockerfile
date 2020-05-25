FROM ubuntu:base as builder1
#this is a temp docker image
LABEL description="Image for leafcutter"
LABEL maintainer="Alexander Dietrich"

ARG star_version=2.7.3a


#install STAR
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean

#get R into image
FROM rocker/tidyverse:3.4.0

#needed to later install devtools
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev

ARG star_version=2.7.3a

#Installing samtools
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y samtools 

#Copying STAR from builder docker
COPY --from=builder1 /docker_main/STAR-${star_version} /opt/STAR-${star_version}/
WORKDIR  /docker_main/

#install stan package, on which leafcutter relies upon
RUN R -e "install.packages('rstan',dependencies=TRUE, repos='http://cran.rstudio.com/')"

#install leafcutter
#RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('davidaknowles/leafcutter/leafcutter')"


ENTRYPOINT ["/bin/bash"]
