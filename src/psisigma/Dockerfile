FROM ubuntu:base
LABEL description="Image for PSI-Sigma"
LABEL maintainer="Tim Faro"

WORKDIR /docker_main

ARG psi_sigma_version=v1.9j
ARG samtools_version=1.11

# Install dependencies
RUN apt-get update && apt-get install --yes build-essential
RUN apt-get install --yes wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libgsl-dev

# Install PSI-Sigma
RUN wget -O - https://github.com/wososa/PSI-Sigma/archive/${psi_sigma_version}.tar.gz | tar -xzv --transform 's:^[^/]*:PSI_Sigma:'

# Install Perl modules
RUN cpan App::cpanminus
RUN cpanm PDL::LiteF
RUN cpanm PDL::Stats
RUN cpanm PDL::GSL::CDF
RUN cpanm Statistics::Multtest

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 && \
    tar xvjf samtools-${samtools_version}.tar.bz2 && \
    cd samtools-${samtools_version} && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    rm -rf samtools-${samtools_version}

RUN mkdir afolder
RUN chmod 777 afolder

#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
