FROM ubuntu:xenial
LABEL description="Image for "
LABEL maintainer="Amit Fenn"

# Needed for smooth OS Install
ARG DEBIAN_FRONTEND=noninteractive

#link to ref genome
# ARG refgenlink=ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

# Binding the folder where the fastq files are. Also needed to restrict the host filesystem within this container in the docker filesystem.
VOLUME /myvol1


# Base OS install.
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y \
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

#Installing tool



COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
