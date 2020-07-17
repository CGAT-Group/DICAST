FROM ubuntu:base as builder1
#this is a temp docker image
LABEL description="Image for JUM"
LABEL maintainer="Alexander Dietrich"

ARG star_version=2.7.3a


#install STAR
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean

#get R into image
FROM rocker/r-base:3.6.3

ARG star_version=2.7.3a

#Installing samtools & bedtools
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y samtools bedtools

#Copying STAR from builder docker
COPY --from=builder1 /docker_main/STAR-${star_version} /opt/STAR-${star_version}/
WORKDIR  /opt/

#install samtools & bedtools
#RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y samtools bedtools

#Perl Modules
RUN cpan install  Array::Utils Statistics::Descriptive

#add binaries to PATH
#ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/STAR-${star_version}/bin/Linux_x86_64



#Install JUM

RUN wget https://github.com/qqwang-berkeley/JUM/releases/download/JUM_v2.0.2/JUM_2.0.2.tar.gz && tar -xzf JUM_2.0.2.tar.gz && rm -rf JUM_2.0.2.tar.gz

#add binaries to PATH
#also add JUM script folder into opt/
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/STAR-${star_version}/bin/Linux_x86_64:/opt/JUM_2.0.2


WORKDIR /docker_main

ENTRYPOINT ["/bin/bash"]

