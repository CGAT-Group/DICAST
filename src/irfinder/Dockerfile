FROM gitlab.lrz.de:5005/ge46ban/dockers/develop/base:0.1 as builder1
#this is a temp docker image
LABEL description="Builder image for STAR"
LABEL maintainer="Alexander Dietrich"

ARG star_version=2.7.4a


#install STAR
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean


#this is the final docker image
FROM gitlab.lrz.de:5005/ge46ban/dockers/develop/base:0.1
LABEL description="Image for IRFinder"
LABEL maintainer="Alexander Dietrich"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

WORKDIR /opt

ARG star_version=2.7.4a
#Copying STAR from builder docker
COPY --from=builder1 /docker_main/STAR-${star_version} /opt/STAR-${star_version}/

#install bedtools
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y bedtools 

#install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && tar -xjvf samtools-1.10.tar.bz2 && cd /opt/samtools-1.10 && ./configure --prefix=/opt/ && make && make install


#install IRFinder from github
RUN wget https://github.com/williamritchie/IRFinder/archive/v1.3.1.tar.gz && tar -zxvf v1.3.1.tar.gz && cd /opt/IRFinder-1.3.1/src/irfinder && make clean && make 
RUN chmod -R 777 IRFinder-1.3.1/bin

#configure PATH variable
ENV PATH=PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/IRFinder-1.3.1/bin/:/opt/samtools-1.10:/opt/STAR-${star_version}/bin/Linux_x86_64

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
