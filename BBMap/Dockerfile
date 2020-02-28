FROM openjdk:11.0.5-stretch
ENV VERSION="38.43"
LABEL description="Image for BBMap "
LABEL maintainer="Amit Fenn"

#link to ref genome
ARG refgenlink=ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Binding the folder where the fastq files are.
VOLUME /myvol1

#Installing tool
RUN cd / && wget "https://sourceforge.net/projects/bbmap/files/BBMap_${VERSION}.tar.gz/download" && \
    tar xf download


#Downloading References
# RUN wget ${refgenlink} -P /docker_main && ref_gen=$(find / -name "*Homo*.fa*") && zcat -d ${ref_gen}>     $(echo "${ref_gen}"|sed s/.gz//)&& rm ${ref_gen}


COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

