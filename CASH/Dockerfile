FROM ubuntu:base
LABEL description="Image for CASH"
LABEL maintainer="Alexander Dietrich"

ARG cash_download=https://sourceforge.net/projects/cash-program/files/2.2.1/cash_v2.2.1.zip

RUN apt-get update && \
	apt-get install -y openjdk-8-jdk && \
  wget ${cash_download} && \
  unzip cash_v2.2.1.zip

WORKDIR /cash_v2.2.1

RUN java -jar cash.jar --help

ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
