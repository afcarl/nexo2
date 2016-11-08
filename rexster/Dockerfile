# Dockerfile for Rexster server 2.4.0 (NeXO)
FROM java:8

MAINTAINER Keiichiro Ono "kono@ucsd.edu"
RUN apt-get update && apt-get install -y wget

RUN mkdir /nexo-rexster
WORKDIR /nexo-rexster

# Download NeXO database
RUN wget http://chianti.ucsd.edu/~kono/ci/data/nexo/nexo-final.tar.gz
# Download rexster v 2.4.0
RUN wget http://chianti.ucsd.edu/~kono/ci/data/atgo/rexster-server-2.4.0.zip

# Setup rexster server
RUN unzip rexster-server-2.4.0.zip

# Add Config file
RUN rm -f /nexo-rexster/rexster-server-2.4.0/config/rexster.xml
ADD ./rexster.xml /nexo-rexster/rexster-server-2.4.0/config/
ADD ./scripts/* /nexo-rexster/rexster-server-2.4.0/scripts/

# Unzip the pre-built database
RUN tar zxvf nexo-final.tar.gz

# Default rexster port
EXPOSE 8182

WORKDIR /nexo-rexster/rexster-server-2.4.0

CMD ["./bin/rexster.sh", "--start"]
# CMD ["/bin/bash"]
