FROM ubuntu:base 
LABEL description="Image for spladder" 
LABEL maintainer="Alexander Dietrich"


#install python3 and pip3
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y python3 python3-dev python3-pip
RUN python3 -m pip install pip --upgrade

#install setuptools & spladder packages 
RUN pip3 install --no-cache-dir -U setuptools
RUN pip3 install --no-cache-dir -U spladder



ENTRYPOINT ["/bin/bash"]
