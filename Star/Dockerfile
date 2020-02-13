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
    && mkdir -p /docker_main/Genomedir /myvol1/Star-Output 
    
WORKDIR /docker_main

# Installing STAR
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean 

ENV PATH=$PATH:/docker_main/STAR-2.7.3a/bin/Linux_x86_64/

#Downloading References
RUN wget ${refgenlink} -P /docker_main && ref_gen=$(find / -name "*Homo*.fa*") && zcat -d ${ref_gen}>     $(echo "${ref_gen}"|sed s/.gz//)&& rm ${ref_gen} \
 && wget ${refantlink} -P /docker_main && ref_ant=$(find / -name "*Homo*.gtf*")&& zcat -d ${ref_ant}>    $(echo "${ref_ant}" |sed s/.gz//)&& rm ${ref_ant} 


#Building the Genome-Index
RUN /docker_main/STAR-${star_version}/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /docker_main/Genomedir --genomeFastaFiles $(find / -name "*Homo*.fa") --runThreadN 120 --sjdbGTFfile $(find / -name "*Homo*.gtf") --sjdbOverhang 100

USER afenn

CMD cd /myvol1 && find /myvol1/ -name "*.fastq"| xargs -n 1 -P 12 /docker_main/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir /docker_main/Genomedir --outFileNamePrefix /myvol1/Star-Output_mapped_ --sjdbGTFfile /docker_main/Homo_sapiens.GRCh38.99.gtf  --twopassMode Basic --runThreadN 40 --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn 



#CMD  ["--genomeDir /docker_main/Genomedir --outFileNamePrefix mapped_ --sjdbGTFfile /docker_main/Homo_sapiens.GRCh38.99.gtf --twopassMode Basic --runThreadN 40 --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn $( find /myvol1/ -name "*.fastq")"]
