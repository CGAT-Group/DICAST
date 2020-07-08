FROM ubuntu:base as builder1
#this is a temp docker image
LABEL description="Builder image for STAR"
LABEL maintainer="Alexander Dietrich"

ARG star_version=2.7.4a


#install STAR
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean


#this is the final docker image
FROM rocker/r-base:4.0.0
LABEL description="Image for KisSplice and KissDE"
LABEL maintainer="Alexander Dietrich"


ARG star_version=2.7.4a
#Copying STAR from builder docker
COPY --from=builder1 /docker_main/STAR-${star_version} /opt/STAR-${star_version}/



#############
# KisSplice #
#############

#install python3
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y python3 python3-dev python3-pip

#download and install KisSplice in /opt
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y cmake zlib1g
WORKDIR /opt
RUN wget ftp://pbil.univ-lyon1.fr/pub/logiciel/kissplice/download/kissplice-2.5.1.tar.gz && tar -xzf /opt/kissplice-2.5.1.tar.gz && rm -f /opt/kissplice-2-5-1.tar.gz &&cd /opt/kissplice-2.5.1 && cmake . && make


########################
#  KisSplice2refgenome #
########################

#kissplice2refgenome requires python2.7 for some reason.......................

#install python 2.7 and pip for that specific version
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y python python-dev && wget https://bootstrap.pypa.io/get-pip.py && python2.7 get-pip.py && rm -f get-pip.py
RUN pip2.7 install biopython==1.76 intervaltree intervaltree-bio
WORKDIR /opt
RUN wget ftp://pbil.univ-lyon1.fr/pub/logiciel/kissplice/tools/kissplice2refgenome-1.2.3.tar.gz && tar -xvf /opt/kissplice2refgenome-1.2.3.tar.gz && rm -f /opt/kissplice2refgenome-1.2.3.tar.gz && cd /opt/kissplice2refgenome-1.2.3 && python2.7 setup.py install --user 

#add Kissplice & Kissplice2refgenome & STAR to path
ENV PATH=PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/kissplice-2.5.1/bin:/opt/kissplice2refgenome-1.2.3:/opt/STAR-${star_version}/bin/Linux_x86_64

#############
#  KissDE   #
#############

WORKDIR /docker_main 
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev 

RUN R -e "install.packages('BiocManager',dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('kissDE',update=TRUE,ask=FALSE)" 



ENTRYPOINT ["/bin/bash"]

