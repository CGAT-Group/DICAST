FROM rocker/r-base:3.6.3
LABEL description="Image for ASpli"
LABEL maintainer="Alexander Dietrich"


RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools


#Installing conda 

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

# Installing Bioconda 
RUN bash -c "source ~/.bashrc && \
conda config --add channels defaults &&\
conda config --add channels bioconda &&\
conda config --add channels conda-forge"

# Installing R packages
RUN R -e "install.packages(c('optparse','BiocManager'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('rtracklayer',update=TRUE,ask=FALSE)"
RUN R -e "BiocManager::install('ASpli',update=TRUE,ask=FALSE)"

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./ASpli.R /docker_main/ASpli.R
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]


