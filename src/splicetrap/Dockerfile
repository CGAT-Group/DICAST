FROM rocker/r-base:3.6.3
LABEL description="Image for SpliceTrap"
LABEL maintainer="Alexander Dietrich"


#install unzip & python (needed for Bowtie I guess...)
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y unzip python-dev samtools

WORKDIR /opt
#install Bowtie
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/bowtie-1.2.3-linux-x86_64.zip && unzip /opt/bowtie-1.2.3-linux-x86_64.zip && rm /opt/bowtie-1.2.3-linux-x86_64.zip && mv /opt/bowtie-1.2.3-linux-x86_64 bowtie-1.2.3

#add bowtie to path
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bowtie-1.2.3


WORKDIR /opt

#download & install splicetrap
RUN wget -O SpliceTrap.0.90.5.tar.gz http://rulai.cshl.edu/splicetrap/download.php?filename=SpliceTrap.0.90.5.tar.gz && tar -xvzf /opt/SpliceTrap.0.90.5.tar.gz && rm /opt/SpliceTrap.0.90.5.tar.gz && cd /opt/SpliceTrap.0.90.5 && echo  | make install 

#add SpliceTrap to path
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/bowtie-1.2.3:/opt/SpliceTrap.0.90.5

WORKDIR /docker_main


ENTRYPOINT ["/bin/bash"]

