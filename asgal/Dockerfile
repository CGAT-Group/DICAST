FROM ubuntu:bioconda
ENV tool=ASGAL
LABEL description="Image for ASGAL"
LABEL maintainer="Amit Fenn"

WORKDIR /docker_main

RUN echo "deb http://in.archive.ubuntu.com/ubuntu/ bionic main restricted universe multiverse \n\
deb http://in.archive.ubuntu.com/ubuntu/ bionic-updates main restricted universe multiverse \n\
deb http://in.archive.ubuntu.com/ubuntu/ bionic-backports main restricted universe multiverse \n\
deb http://security.ubuntu.com/ubuntu bionic-security main restricted universe multiverse \n\
deb-src http://in.archive.ubuntu.com/ubuntu/ bionic main restricted universe multiverse \n\
deb-src http://security.ubuntu.com/ubuntu bionic-security main restricted universe multiverse \n\
deb-src http://in.archive.ubuntu.com/ubuntu/ bionic-backports main restricted universe multiverse \n\
deb-src http://in.archive.ubuntu.com/ubuntu/ bionic-updates main restricted universe multiverse " >>  /etc/apt/sources.list


#Installing dependencies 
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y \
	python3 \
	python3-pip \
	python3-setuptools \
	python3-biopython \
	python3-biopython-sql \
	python3-pysam \
	cmake \
	libboost1.65-all-dev \
	samtools \
	unzip \
	wget \
	curl \
	zlib1g-dev \
	liblzma-dev \
	libjemalloc-dev \
	libjemalloc1 \
	libghc-bzlib-dev \
	libgff-dev \
	libtbb-dev; \
    pip3 install gffutils

#Installing ASGAL
RUN cd /opt/ ;\
git clone --recursive https://github.com/AlgoLab/galig.git ;\
cd galig/salmon ;\
wget https://github.com/AlgoLab/galig/files/4437983/CMakeLists.txt.patch.txt ;\
git apply CMakeLists.txt.patch.txt ;\
cd /opt/galig/ ;\
make salmon ;\
make lemon sdsl ;\
make

RUN conda install -c anaconda biopython -y

ENV PATH=$PATH:/opt/galig:/opt/galig/salmon/bin

#ENTRYPOINT SCRIPT
# COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
#RUN chmod 777 /docker_main/ENTRYPOINT.sh
# ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
CMD ["/bin/bash"]
