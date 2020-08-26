FROM ubuntu:base
LABEL description="Image for MAJIQ"
LABEL maintainer="Alexander Dietrich"

WORKDIR /opt
#install htslib
RUN cd /usr/bin
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar -vxjf htslib-1.9.tar.bz2 && rm -f htslib-1.9.tar.bz2 && cd htslib-1.9 && make && make install

#install python3.6 and update pip
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa -y -u
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y python3.6 python3.6-dev python3-pip python3.6-venv
RUN python3.6 -m pip install pip --upgrade

#install majiq dependencies
RUN pip3.6 install --no-cache-dir -U wheel setuptools 
RUN pip3.6 install --no-cache-dir -U numpy cython
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y git

#download and install majiq from source
WORKDIR /tmp
RUN git clone https://bitbucket.org/biociphers/majiq_stable.git
WORKDIR /tmp/majiq_stable
RUN pip3.6 install --no-cache-dir -r requirements.txt
RUN python3.6 setup.py install
RUN rm -rf /tmp
WORKDIR /opt

#download and add script to convert GTF file to GFF3 format
RUN wget --no-check-certificate http://majiq.biociphers.org/download/gtf2gff3.pl && wget http://www.sequenceontology.org/resources/gtf2gff3.cfg

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt:/opt/htslib-1.9

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

