FROM bowtie:0.2

LABEL software="ContextMap"
LABEL maintainer="Amit Fenn"

# Name of mounted folder
VOLUME /MOUNT

#Copying Java runtime install file
#COPY [--chown=nobody:nogroup] ./jre-8u241-linux-i586.tar.gz /home/biodocker/jre-8u241-linux-i586.tar.gz

WORKDIR /home/biodocker/
RUN wget https://javadl.oracle.com/webapps/download/AutoDL?BundleId=242978_a4634525489241b9a9e1aa73d9e118e6 && mv AutoDL\?BundleId\=242978_a4634525489241b9a9e1aa73d9e118e6 jre-8u241-linux-i586.tar.gz

#Installations
RUN wget https://www.bio.ifi.lmu.de/software/contextmap/contextmap_source_v2_7_9.zip --no-check-certificate && unzip contextmap_source_v2_7_9.zip && tar -zxvf jre-8u241-linux-i586.tar.gz && rm contextmap_source_v2_7_9.zip jre-8u241-linux-i586.tar.gz


COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /home/biodocker/ENTRYPOINT.sh
#RUN chmod +x /home/biodocker/ENTRYPOINT.sh
#ENTRYPOINT ["/home/biodocker/ENTRYPOINT.sh"]

CMD ["/bin/bash","/home/biodocker/ENTRYPOINT.sh"]
