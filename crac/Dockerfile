FROM ubuntu:base
LABEL description="Image for "
LABEL maintainer="Amit Fenn"

# Needed for smooth OS Install
ARG DEBIAN_FRONTEND=noninteractive

VOLUME /myvol1
    
WORKDIR /docker_main

#installing dependencies
RUN apt-get --no-install-recommends --fix-broken install -y autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

#Installing tool
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && wget https://gforge.inria.fr/frs/download.php/file/35190/crac-2.4.0.tar.gz && unar htslib-1.10.2.tar.bz2 && tar -xvzf crac-2.4.0.tar.gz && cd /docker_main/htslib-1.10.2/ && ./configure --prefix=/usr/ && make && make install && cd /docker_main/crac-2.4.0/ && ./configure && make && make install 
 


COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
