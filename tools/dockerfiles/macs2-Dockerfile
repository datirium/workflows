#################################################################
# Dockerfile
#
# Software:         MACS2
# Software Version: 2.1.1.20160309
# Description:      MACS2 image for BioWardrobe
# Website:          http://liulab.dfci.harvard.edu/MACS/
# Provides:         MACS2
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/macs2:v2.1.1 -f macs2-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/macs2:v2.1.1
# Run Cmd:          docker run --rm -ti biowardrobe2/macs2:v2.1.1 macs2 -h
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing MACS2

ENV VERSION 2.1.1.20160309
ENV NAME MACS2

RUN pip install numpy && \
    pip install wheel && \
    pip install setuptools && \
    pip install ${NAME}==${VERSION}
