FROM dicast/base:0.1
ENV tool=osa
LABEL description="Image for $tool"
LABEL maintainer="Amit Fenn"

WORKDIR /opt

#Installing dependencies
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF && apt install apt-transport-https ca-certificates -y && echo "deb https://download.mono-project.com/repo/ubuntu stable-xenial main" > /etc/apt/sources.list.d/mono-official-stable.list
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y \
libgdiplus zlib1g-dev mono-devel 
RUN sed -i '46i <dllmap dll="gdiplus.dll" target="/usr/lib/libgdiplus.so.0"/>' /etc/mono/config 

#Installing OSHELL
RUN mkdir ~/oshell && cd ~/oshell && wget -c  http://omicsoft.com/software_update/OmicsoftUpdater.exe && touch oshell.exe && /usr/bin/mono OmicsoftUpdater.exe






#ENTRYPOINT SCRIPT
#COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
#RUN chmod 777 /docker_main/ENTRYPOINT.sh
#ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
