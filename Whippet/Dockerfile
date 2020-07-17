FROM ubuntu:base
LABEL description="Image for Whippet"
LABEL maintainer="Alexander Dietrich"


WORKDIR /opt
#download and install julia; need to change folder name manually, since it is not consistent
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz && tar -xvzf julia-0.6.4-linux-x86_64.tar.gz && rm julia-0.6.4-linux-x86_64.tar.gz && mv julia-* julia-0.6.4

#export julia binary to path
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/julia-0.6.4/bin
ENV JULIA_PKGDIR=/opt/.julia

#install whippet and add to path
RUN julia -e 'Pkg.init(); Pkg.resolve(); Pkg.add("Whippet"); using Whippet'
#RUN cp ~/.julia/v0.6/REQUIRE /opt/.julia/v0.6 && julia -e 'Pkg.resolve()'
#RUN julia -e 'Pkg.add("Whippet"); using Whippet' 
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/julia-0.6.4/bin:/opt/.julia/v0.6/Whippet

WORKDIR /docker_main
RUN ln -s /opt/.julia/v0.6/Whippet/bin /docker_main
RUN chmod 777 /opt/.julia/lib/v0.6/*


ENTRYPOINT ["/bin/bash"]
