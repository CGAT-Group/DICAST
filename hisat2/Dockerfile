FROM ubuntu:base
LABEL description="Image for hisat2" 
LABEL maintainer="Amit Fenn"

WORKDIR /docker_main

# Installing hisat2 
RUN wget http://ccb.jhu.edu/software/hisat2/downloads/hisat2-2.0.0-beta-source.zip  && unar hisat2-2.0.0-beta-source.zip && cd /docker_main/hisat2-2.0.0-beta/ && make 

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
