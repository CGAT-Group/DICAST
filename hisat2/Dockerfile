FROM ubuntu:xenial
LABEL description="Image for hisat2" 
LABEL maintainer="Amit Fenn"

# Needed for smooth OS Install
ARG DEBIAN_FRONTEND=noninteractive

#link to ref genome
ARG refgenlink=ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
#link to ref annotation
ARG refantlink=ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz 

# Binding the folder where the fastq files are.
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
    zlib1g-dev\
    unar \
    && mkdir -p /myvol1/hisat-Output 
    
WORKDIR /docker_main

# Installing hisat2 
RUN wget http://ccb.jhu.edu/software/hisat2/downloads/hisat2-2.0.0-beta-source.zip  && unar hisat2-2.0.0-beta-source.zip && cd /docker_main/hisat2-2.0.0-beta/ && make 

ENV PATH=$PATH:/docker_main/hisat2-2.0.0-beta/

#Downloading References
RUN wget ${refgenlink} -P /docker_main && ref_gen=$(find / -name "*Homo*.fa*") && zcat -d ${ref_gen}>     $(echo "${ref_gen}"|sed s/.gz//)&& rm ${ref_gen} 
# && wget ${refantlink} -P /docker_main && ref_ant=$(find / -name "*Homo*.gtf*")&& zcat -d ${ref_ant}>    $(echo "${ref_ant}" |sed s/.gz//)&& rm ${ref_ant} 

# Building genome index	
# RUN hisat2-build /docker_main/Homo_sapiens.GRCh38.dna.primary_assembly.fa latest_human && python /docker_main/hisat2-2.0.0-beta/extract_splice_sites.py /docker_main/Homo_sapiens.GRCh38.99.gtf > splicesites.txt


COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
