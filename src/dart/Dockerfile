FROM gitlab.lrz.de:5005/ge46ban/dockers/develop/base:0.1
ENV tool=dart
LABEL description="Image for dart "
LABEL maintainer="Amit Fenn"

WORKDIR /opt

#Installing dependencies
RUN apt-get update -y && apt-get install libboost-all-dev libbz2-dev liblzma-dev --no-install-recommends --fix-broken -y 

#Installing tool
#RUN git clone https://github.com/hsinnan75/Dart.git v1.4.0 ;\
RUN git clone https://github.com/hsinnan75/DART ;\
cd /opt/DART ;\
make

ENV PATH=$PATH:/opt/DART/bin

# Name of mounted folder
VOLUME /MOUNT


#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
