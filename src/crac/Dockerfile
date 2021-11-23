FROM gitlab.lrz.de:5005/ge46ban/dockers/develop/crac:0.3

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]