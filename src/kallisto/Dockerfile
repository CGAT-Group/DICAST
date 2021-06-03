FROM ubuntu:base
ENV tool=kallisto
LABEL description="Image for $tool"
LABEL maintainer="Amit Fenn"

WORKDIR /docker_main

#Installing tool
WORKDIR /opt
RUN git clone https://github.com/pachterlab/kallisto.git
WORKDIR kallisto

#Installing htslib
RUN git clone https://github.com/samtools/htslib
RUN rm -rf -f build
RUN rm -rf /ext/htslib
RUN cp -r htslib /ext/
WORKDIR /opt/kallisto/ext/htslib
RUN autoconf && autoheader

#Kallisto
WORKDIR /opt/kallisto
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make
RUN make install

COPY Dockerfile /opt/.

#ENTRYPOINT SCRIPT
# COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
# RUN chmod 777 /docker_main/ENTRYPOINT.sh
# ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
