FROM bowtie:0.1

LABEL software="ContextMap"
LABEL maintainer="Amit Fenn"

#Copying Java runtime install file
COPY [--chown=nobody:nogroup] ./jre-8u241-linux-i586.tar.gz /home/biodocker/jre-8u241-linux-i586.tar.gz

WORKDIR /home/biodocker/

#Installations
RUN wget https://www.bio.ifi.lmu.de/software/contextmap/contextmap_source_v2_7_9.zip && unzip contextmap_source_v2_7_9.zip && tar -zxvf jre-8u241-linux-i586.tar.gz && rm contextmap_source_v2_7_9.zip jre-8u241-linux-i586.tar.gz

USER root
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /home/biodocker/ENTRYPOINT.sh
RUN chmod +x /home/biodocker/ENTRYPOINT.sh
#ENTRYPOINT ["/home/biodocker/ENTRYPOINT.sh"]

CMD ["/bin/bash","/home/biodocker/ENTRYPOINT.sh"]
