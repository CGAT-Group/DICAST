FROM ubuntu:base
ENV tool=dart
LABEL description="Image for dart "
LABEL maintainer="Amit Fenn"

WORKDIR /opt

#Installing dependencies
# RUN apt-get update -y && apt-get --no-install-recommends --fix-broken inst    all -y \

#Installing tool
RUN git clone https://github.com/hsinnan75/Dart.git v1.4.0 ;\
cd /opt/v1.4.0/ ;\
make

ENV PATH=$PATH:/opt/v1.4.0/bin


#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
