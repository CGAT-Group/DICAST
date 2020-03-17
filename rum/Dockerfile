FROM ubuntu:base 
ENV tool=rum	
LABEL description="Image for $tool"
LABEL maintainer="Amit Fenn"

WORKDIR /opt/

#Installing tool
RUN wget https://github.com/PGFI/rum/archive/master.tar.gz && tar -zxvf master.tar.gz && rm master.tar.gz && cd rum-master/ &&  perl Makefile.PL && make && make install

WORKDIR /docker_main

#ENTRYPOINT SCRIPT
# COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
# RUN chmod 777 /docker_main/ENTRYPOINT.sh
# ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
