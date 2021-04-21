FROM dicast/base:conda0.1
LABEL description="Image for rMATS"
LABEL maintainer="Tim Faro"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV LD_LIBRARY_PATH /opt/conda/lib/:$LD_LIBRARY_PATH

WORKDIR /docker_main

# Installing tool
RUN conda install -c bioconda rmats -y

#Create temporary folder
RUN mkdir tmp
RUN chmod 777 tmp

#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
