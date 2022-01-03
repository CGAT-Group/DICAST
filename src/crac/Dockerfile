FROM dicastproj/dicast:crac

COPY [--chown=nobody:nogroup] ./ENTRYPOINT.sh /docker_main/
RUN chmod 777 /docker_main/ENTRYPOINT.sh
ENTRYPOINT ["/docker_main/ENTRYPOINT.sh"]