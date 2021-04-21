FROM gitlab.lrz.de:5005/ge46ban/dockers/develop/base-conda:0.1
LABEL description="Image for minimap2" 
LABEL maintainer="Amit Fenn"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# Name of mounted folder
VOLUME /MOUNT

WORKDIR /docker_main

# RUN conda install 
RUN conda install -y minimap2 

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
