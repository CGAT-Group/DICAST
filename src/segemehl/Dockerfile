FROM ubuntu:base
ENV tool=segemehl
LABEL description="Image for segemehl"
LABEL maintainer="Amit Fenn"

WORKDIR /docker_main

#installing htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && unar htslib-1.10.2.tar.bz2 && cd /docker_main/htslib-1.10.2/ &&  ./configure --prefix=/usr/ && make && make install && rm /docker_main/htslib*.bz2
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y \
 pkg-config

WORKDIR /docker_main

#Installing tool
ENV SEGEMEHL_VERSION="0.3.4"
ENV SEGEMEHL_HOME="/opt/segemehl-0.3.4"
ENV PATH="$SEGEMEHL_HOME:$PATH"

RUN wget https://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/segemehl-${SEGEMEHL_VERSION}.tar.gz && tar -xzvf ./s*.gz && rm ./s*.gz && \
    pwd

RUN cd /docker_main/segemehl* && \
    make all && \
    mv /docker_main/segemehl-${SEGEMEHL_VERSION}/ /opt/ && \
    rm -rf /docker_main/segemehl-${SEGEMEHL_VERSION} && \
    apt-get purge -y --auto-remove curl &&  \
    apt-get clean -yq

# Name of mounted folder
VOLUME /MOUNT  

#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
