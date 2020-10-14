FROM ubuntu:base
LABEL description="Image for dSpliceType"
LABEL maintainer="Fanny Roessler"

ENV r_version 3.6.0

RUN apt-get update -y && apt-get install -y \
  gfortran \
  libreadline-dev \
  libpcre3-dev \
  libcurl4-openssl-dev \
  build-essential \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  openjdk-8-jdk \
  wget \
  libssl-dev \
  libxml2-dev \
  libnss-sss \
  git \
  build-essential \
  cmake \
  python3

################################################################################
####################### Install R ##############################################

# change working dir
WORKDIR /usr/local/bin

# install R
RUN wget https://cran.r-project.org/src/base/R-3/R-${r_version}.tar.gz
RUN tar -zxvf R-${r_version}.tar.gz
WORKDIR /usr/local/bin/R-${r_version}
RUN ./configure --prefix=/usr/local/ --with-x=no
RUN make
RUN make install

# install R packages
RUN R --vanilla -e 'install.packages(c("data.table", "plyr", "tidyverse"), repos = "http://cran.us.r-project.org")'


#link to zip folder
ARG dSpliceType_download=https://sourceforge.net/projects/dsplicetype/files/dSpliceType-2.0.0.zip
ARG bedtools_download=https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
#ARG bedtools_download=https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary
ARG regtools_download=https://github.com/griffithlab/regtools.git

# Name of mounted volume
VOLUME /myvol1

WORKDIR /opt

################################################################################
##################### Install Regtools #########################################

# clone git repository
RUN cd / && git clone https://github.com/griffithlab/regtools.git

# make a build directory for regtools
WORKDIR /regtools/


# compile from source
RUN mkdir build && cd build && cmake .. && make


#RUN apt-get update && wget ${dSpliceType_download} && unzip dSpliceType-2.0.0.zip && rm dSpliceType-2.0.0.zip && chmod +x dSpliceType.jar
#RUN apt-get update && wget ${bedtools_download} && tar -zxvf bedtools-2.29.1.tar.gz && cd bedtools2 && make && cp bin/* /usr/local/bin/
#RUN apt-get update && wget ${bedtools_download} && mv bedtools.static.binary bedtools && chmod a+x bedtools && bedtools coverage
RUN apt-get update && git clone ${regtools_download} && cd regtools/ && mkdir build && cd build/ && cmake .. && make

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /$HOME


################################################################################
###################### set environment path    #################################

# make a build directory for regtools
WORKDIR /regtools/scripts/

# add regtools executable to path
ENV PATH="/regtools/build:/usr/local/bin/R-${r_version}:${PATH}"


#COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
#RUN chmod 777 /docker_main/ENTRYPOINT.sh
#ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

