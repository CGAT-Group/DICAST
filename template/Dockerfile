FROM ubuntu:base	#build in this git repo under BaseOS
ENV tool=		# insert name of tool-here
LABEL description="Image for $tool"
LABEL maintainer="Amit Fenn"

WORKDIR /docker_main

#Installing dependencies
# RUN apt-get update -y && apt-get --no-install-recommends --fix-broken inst    all -y \

#Installing tool
RUN 

#ENTRYPOINT SCRIPT
COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/ENTRYPOINT.sh
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]
