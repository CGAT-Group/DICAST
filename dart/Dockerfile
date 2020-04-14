FROM ubuntu:base
ENV tool=dart
LABEL description="Image for $tool"
LABEL maintainer="Amit Fenn"

WORKDIR /opt

#Installing dependencies
# RUN apt-get update -y && apt-get --no-install-recommends --fix-broken inst    all -y \

#Installing tool
RUN git clone https://github.com/hsinnan75/Dart.git ;\
cd /opt/Dart/ ;\
git fetch ;\
git merge origin/master master ;\
sed -i '823,824 s/^/\/\//' /opt/Dart/src/Mapping.cpp ;\
make

ENV PATH=$PATH:/opt/Dart/bin



#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
