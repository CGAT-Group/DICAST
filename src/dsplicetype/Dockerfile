FROM ubuntu:base
LABEL description="Image for dSpliceType"
LABEL maintainer="Fanny Roessler"

ARG dSpliceType_download=https://sourceforge.net/projects/dsplicetype/files/dSpliceType-2.0.0.zip
ARG regtools_download=https://github.com/griffithlab/regtools.git

### install dSpliceType ###
RUN apt-get update && wget ${dSpliceType_download} && unzip dSpliceType-2.0.0.zip && rm dSpliceType-2.0.0.zip && chmod +x dSpliceType.jar


### install regtools ###
# clone git repository
RUN cd / && git clone ${regtools_download}

# make a build directory for regtools
WORKDIR /regtools/

# compile from source
RUN mkdir build && cd build && cmake .. && make


### install bedtools ###
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y bedtools 


################################## docker setup ###

ENV PATH="/regtools/build:/usr/bin/bedtools:${PATH}"

WORKDIR /docker_main/

# Name of mounted folder
VOLUME /MOUNT


COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

