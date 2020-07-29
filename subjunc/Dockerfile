FROM ubuntu:base
ENV tool=subjunc
LABEL description="Image for subjunc"
LABEL maintainer="Amit Fenn"

WORKDIR /opt/

#Installing tool
RUN wget wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz/download -O subread-2.0.0-Linux-x86_64.tar.gz ; tar -xvzf subread-2.0.0-Linux-x86_64.tar.gz 
# && rm /opt/subread-2.0.0-Linux-x86_64.tar.gz

RUN ls /opt/subread-2.0.0-Linux-x86_64

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

