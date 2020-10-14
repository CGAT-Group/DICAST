FROM openjdk:11.0.5-stretch
ENV VERSION="38.43"
LABEL description="Image for BBMap "
LABEL maintainer="Amit Fenn"

#link to ref genome
ARG refgenlink=ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Name of mounted folder
VOLUME /MOUNT

#Installing tool
RUN cd / && wget "https://sourceforge.net/projects/bbmap/files/BBMap_${VERSION}.tar.gz/download" && \
    tar xf download

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
