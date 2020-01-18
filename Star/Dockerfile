FROM ubuntu:xenial

LABEL description="Image for STAR aligner version 2.7.3a"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y --no-install-recommends \
build-essential \
    bzip2 \
    cmake \
    default-jdk \
    git \
    libnss-sss \
    libtbb2 \
    libtbb-dev \
    ncurses-dev \
    python-dev \
    python-pip \
    tzdata \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev

WORKDIR /docker_main
RUN wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz --no-check-certificate && tar -xzf 2.7.3a.tar.gz && cd STAR-2.7.3a/source && make STAR


RUN cd /docker_main / && \
   rm -rf ${star_version}.tar.gz && \
   apt-get autoremove -y && \
   apt-get autoclean -y  && \
   apt-get clean

# Set default working path
WORKDIR /docker_main

ENTRYPOINT {"bash"}
