FROM ubuntu:base
LABEL description="Image for Whippet"
LABEL maintainer="Alexander Dietrich"

#install samtools
RUN apt-get update && apt-get --no-install-recommends --fix-broken -y install samtools


WORKDIR /opt
#download and install julia; need to change folder name manually, since it is not consistent
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz && tar -xvzf julia-0.6.4-linux-x86_64.tar.gz && rm julia-0.6.4-linux-x86_64.tar.gz && mv julia-* julia-0.6.4

#export julia binary to path
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/julia-0.6.4/bin
#install julia packages into /opt
ENV JULIA_PKGDIR=/opt/.julia
ENV JULIA_HISTORY=/opt/.julia_history

#install whippet 
RUN julia -e 'Pkg.init(); Pkg.resolve(); Pkg.add("Whippet"); using Whippet'

WORKDIR /docker_main
#create link to whippet scripts
RUN ln -s /opt/.julia/v0.6/Whippet/bin /docker_main/
#run script once in "dry" mode to precompile all libraries and change permissions of them
RUN julia /docker_main/bin/whippet-index.jl -h
RUN julia /docker_main/bin/whippet-delta.jl -h
RUN chmod 777 -R /opt/.julia/lib/v0.6/*


COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]

