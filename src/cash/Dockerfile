FROM gitlab.lrz.de:5005/ge46ban/dockers/develop/base:0.1
LABEL description="Image for CASH"
LABEL maintainer="Alexander Dietrich"

#install samtools
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y samtools


#link to zip folder
ARG cash_download=https://sourceforge.net/projects/cash-program/files/2.2.1/cash_v2.2.1.zip

WORKDIR /opt

#download and install cash
RUN apt-get update && wget ${cash_download} && unzip cash_v2.2.1.zip && rm cash_v2.2.1.zip && chmod +x cash_v2.2.1/cash.jar

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/cash_v2.2.1

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

