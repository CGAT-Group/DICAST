FROM ubuntu:base
LABEL description="Image for crac"
LABEL maintainer="Amit Fenn"
   
WORKDIR /docker_main

#Installing dependencies: htslib.
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && unar htslib-1.10.2.tar.bz2 && cd /docker_main/htslib-1.10.2/ &&  ./configure --prefix=/usr/ && make && make install

# Name of mounted folder
VOLUME /MOUNT

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

#Installing tool
#	RUN wget https://gforge.inria.fr/frs/download.php/file/35191/crac_2.4.0_amd64.deb -p /docker_main/ && dpkg -i /docker_main/crac_2.4.0_amd64.deb && apt-get install -f
ARG crac_version="2.5.2"
RUN wget https://gforge.inria.fr/frs/download.php/latestfile/3502/crac-${crac_version}.tar.gz --no-check-certificate && unar crac-${crac_version}.tar.gz && cd crac-${crac_version} && ./configure && make && echo check && make check  && make install


COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
