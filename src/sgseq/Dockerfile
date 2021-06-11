FROM  dicast/star:0.2 as builder

FROM rocker/r-base:4.0.0
LABEL description="Image for SGSeq"
LABEL maintainer="Alexander Dietrich"

ARG star_version=2.7.5c
#Installing conda 


#Copying STAR from builder docker
COPY --from=builder /docker_main/STAR-${star_version} /opt/STAR-${star_version}/
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/STAR-${star_version}/bin/Linux_x86_64
RUN apt update && apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev

# Installing conda
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
conda config --add channels defaults && \
conda config --add channels bioconda && \
conda config --add channels conda-forge"


# Installing SGSeq
RUN R -e "install.packages(c('BiocManager','optparse'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('SGSeq','GenomicFeatures'),update=TRUE,ask=FALSE)" 
RUN chmod 777 /usr/local/ -R 

RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y \
 python3.8 python3-pip

WORKDIR /docker_main

COPY [--chown=nobody:nogroup] ./SGSeq.R /docker_main/SGSeq.R
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]