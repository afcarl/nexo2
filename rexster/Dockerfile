# Dockerfile for Neo4j

# Neo4j Server
# Repository http://github.com/neo4j-contrib/docker-neo4j
FROM java:openjdk-8-jre

MAINTAINER Keiichiro Ono "kono@ucsd.edu"

# Install unzip so we later can unpack
RUN apt-get -y install unzip

WORKDIR /

# Copy last version of Rexter
ADD http://tinkerpop.com/downloads/rexster/rexster-server-2.6.0.zip /

# Add config file

RUN unzip rexster-server-2.6.0.zip
RUN rm rexster-server-2.6.0.zip
RUN mv rexster-server-2.6.0 rexster

RUN mkdir /rexster/scripts
ADD ./scripts /rexster/scripts/

ADD ./rexster.xml /rexster/config/

# Add database
RUN mkdir /webapp && mkdir /webapp/database && mkdir /webapp/database/nexo-dag
ADD ./nexo-neo4j-db /webapp/database/nexo-dag/

# Server listening port
EXPOSE 8182 8184 8183

WORKDIR /rexster

## Entrypoint
CMD ["bin/rexster.sh", "--start"]
