FROM ubuntu:base
LABEL description="Image for CASH"
LABEL maintainer="Alexander Dietrich"

#link to zip folder
ARG cash_download=https://sourceforge.net/projects/cash-program/files/2.2.1/cash_v2.2.1.zip

RUN apt-get update && wget ${cash_download} && unzip cash_v2.2.1.zip

WORKDIR /docker_main/cash_v2.2.1


ENTRYPOINT ["/bin/bash"]

