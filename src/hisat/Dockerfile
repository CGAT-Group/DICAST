FROM dicast/base:conda0.1
LABEL description="Image for hisat2"
LABEL maintainer="Amit Fenn"

WORKDIR /docker_main

# Installing hisat2 v2.2.1
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2.2.2.1.zip && unzip hisat2.2.2.1.zip && cd /docker_main/hisat2-2.2.1/ && ln -s /docker_main/hisat2-2.2.1/* /usr/local/bin/

# Name of mounted folder
VOLUME /MOUNT

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
