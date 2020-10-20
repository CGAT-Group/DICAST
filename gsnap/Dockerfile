FROM ubuntu:base
ENV tool=gsnap
LABEL description="Image for gsnap"
LABEL maintainer="Amit Fenn"

WORKDIR /docker_main

#Installing tool
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-03-12.tar.gz && tar -xvf /docker_main/gmap-gsnap-2020-03-12.tar.gz && rm /docker_main/gmap-gsnap-2020-03-12.tar.gz  && cd gmap-2020-03-12 && mkdir /docker_main/gsnap/gmapdb -p
RUN ./g*/configure --with-gmapdb=/docker_main/gsnap/gmapdb && make && make install

# Name of mounted folder
VOLUME /MOUNT

#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
