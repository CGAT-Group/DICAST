FROM gitlab.lrz.de:5005/ge46ban/dockers/develop/base:0.1
LABEL description="Image for juncBASE"
LABEL maintainer="Alexander Dietrich"

WORKDIR /docker_main

#installs R 3.2.3
RUN apt-get update && apt-get install -y r-base

#install all required python packages (with their versions working with python2.7)
RUN python -m pip install pip --upgrade
RUN pip install --no-cache-dir -U wheel setuptools numpy
RUN pip install --no-cache-dir -U biopython==1.59 pysam rpy2==2.8.6

#install sqlite
RUN apt-get update && apt-get install -y sqlite3 libsqlite3-dev
RUN pip install --no-cache-dir -U pysqlite


#download juncBASE from github
RUN git clone https://github.com/anbrooks/juncBASE.git
WORKDIR /docker_main/juncBASE/

ENTRYPOINT ["/bin/bash"]

