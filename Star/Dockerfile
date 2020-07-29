FROM ubuntu:base as builder
LABEL description="Image for STAR aligner version 2.7.3a"
LABEL maintainer="Amit Fenn"

ARG star_version=2.7.3a
RUN wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz --no-check-certificate -P /docker_main && tar -xzf /docker_main/${star_version}.tar.gz && cd STAR-${star_version}/source && make STAR && rm -rf /docker_main/${star_version}.tar.gz && apt autoremove && apt clean 

#Multistage to decrease size
FROM alpine:latest
ARG star_version=2.7.3a

#Installing Bash
RUN apk update && apk upgrade && apk add bash && apk add bash-completion && echo Bash installed
COPY --from=builder /docker_main/STAR-${star_version} /opt/STAR-${star_version}/



# Installing STAR

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
