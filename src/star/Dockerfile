FROM dicast/base:0.1 as builder
LABEL description="Image for STAR aligner version 2.7.5c"
LABEL maintainer="Amit Fenn"

ARG star_version=2.7.5c
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean 

# Name of mounted folder
VOLUME /MOUNT

# Installing STAR

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

