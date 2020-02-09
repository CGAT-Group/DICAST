FROM ubuntu:xenial
LABEL description="Image for STAR aligner version 2.7.3a"
LABEL maintainer="Amit Fenn"


# Needed for smooth OS Install
ARG DEBIAN_FRONTEND=noninteractive

# Star Version
ARG star_version=2.7.3a

#link to ref genome
ARG refgenlink=ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
#link to ref annotation
ARG refantlink=ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz 

# Binding the folder where the fastq files are.
VOLUME /myvol1

RUN mkdir /docker_main /docker_main/Genomedir /myvol1/Star-Output
WORKDIR /docker_main

# Base OS install.
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

# Installing STAR
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR

#Downloading References
RUN wget ${refgenlink} -P /docker_main && ref_gen=$(find / -name "*Homo*.fa*") && zcat -d ${ref_gen}>     $(echo "${ref_gen}"|sed s/.gz//)&& rm ${ref_gen}
RUN wget ${refantlink} -P /docker_main && ref_ant=$(find / -name "*Homo*.gtf*")&& zcat -d ${ref_ant}>    $(echo "${ref_ant}" |sed s/.gz//)&& rm ${ref_ant} 


#Building the Genome-Index
RUN /docker_main/STAR-${star_version}/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /docker_main/Genomedir --genomeFastaFiles $(find / -name "*Homo*.fa") --runThreadN 4 --sjdbGTFfile $(find / -name "*Homo*.gtf") --sjdbOverhang 100

#Cleaning up
RUN rm -rf /docker_main/${star_version}.tar.gz && \
    apt-get autoremove -y && \
    apt-get autoclean -y  && \
    apt-get clean

#ENTRYPOINT ["/bin/bash"]
ENTRYPOINT ["/docker_main/Star-${star_version}/bin/Linux_x86_64/STAR --genomeDir /docker_main/Genomedir --outFileNamePrefix mapped_ --sjdbGTFfile $(find / -name "*Homo*.gtf*") --twopassMode Basic --runThreadN 40  --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS "]
CMD ["--readFilesIn $(find /myvol1 -name "*.fastq""]
