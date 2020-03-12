FROM ubuntu:base	
LABEL description="Image for dart"
LABEL maintainer="Amit Fenn"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

VOLUME /myvol1

WORKDIR /docker_main

#Installing tool

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

# Create the environment:
COPY dartcondainstall.yml .

# Make RUN commands use the new environment:
# SHELL ["conda", "run", "-n", "dart", "/bin/bash", "-c"]

# RUN conda init bash
RUN conda env create -f dartcondainstall.yml
# RUN conda activate dart
RUN conda env list
RUN conda list

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
# ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
