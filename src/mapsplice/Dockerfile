FROM ubuntu:conda
LABEL description="Image for mapsplice"
LABEL maintainer="Amit Fenn"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# Name of mounted folder
VOLUME /MOUNT

WORKDIR /docker_main

# RUN conda install 
RUN conda install mapsplice -y

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
