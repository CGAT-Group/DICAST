FROM ubuntu:xenial
LABEL description="Image for project base OS"
LABEL maintainer="Amit Fenn"

# Needed for smooth OS Install
ARG DEBIAN_FRONTEND=noninteractive

#link to ref genome
ARG refgenlink=ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

# Binding the folder where the fastq files are.
VOLUME ["/myvol1", "/myvol2"]


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
    zlib1g-dev\
    libbz2-dev \
    liblzma-dev\
    tmux\
    vim\
    unar\
    gcc-multilib \
    apt-utils \
    zlib1g-dev \
    libhdf5-dev \
    git-all \
    autoconf \
    automake \
    libcurl4-openssl-dev
    
#Adding a signature
RUN echo '[ ! -z "$TERM" -a -r /etc/motd ] && cat /etc/issue && cat /etc/motd' \
    >> /etc/bash.bashrc \
    ; echo "\
=========================================================================================\n\
=  Welcome to a Docker container built by the Alternative Splicing - Docker Team at TUM =\n\
=========================================================================================\n\
\n\
If you need help with this Docker, please contact:  ^_^ \n\
Amit FENN, Chair of Experimental Bioinformatics, Technical University of Munich (TUM), Freising, Germany \n\
\n\
The Docker team: Fanny Rößler, Alexander Dietrich and Amit\
\n\
If you liked our Docker, then consider buying us a beer.\n"\
    > /etc/motd

    
WORKDIR /docker_main

#Downloading References
# RUN wget ${refgenlink} -P /docker_main && ref_gen=$(find / -name "*Homo*.fa*") && zcat -d ${ref_gen}>     $(echo "${ref_gen}"|sed s/.gz//)&& rm ${ref_gen} 
 
ENTRYPOINT ["/bin/bash"]
