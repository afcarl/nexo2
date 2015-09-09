#
# Dockerfile for NeXO Term Enrichment service.
#

FROM ipython/scipystack

MAINTAINER Keiichiro Ono <kono@ucsd.edu>

# Install git client
RUN apt-get update && apt-get install -y git

# Required Python dependencies
RUN pip2 install flask-restful

# Form a set of standard directories.
RUN mkdir -p /downloads
#RUN mkdir -p /work

# Install python libraries

# qvalue library
WORKDIR /downloads

# Install qvalue library 
RUN git clone https://github.com/idekerlab/qvalue.git
WORKDIR /downloads/qvalue
RUN python setup.py install

# Add actual application to the container
ADD . /downloads/nexo-term-enrichment
WORKDIR  /downloads/nexo-term-enrichment

# EXPOSE PORT
EXPOSE 5000

RUN python -c "from flask import Flask"
CMD ["python", "./nexo.py"]