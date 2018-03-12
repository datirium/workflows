#################################################################
# Dockerfile
#
# Software:         fastqc
# Software Version: 0.11.5
# Description:      fastqc image for BioWardrobe
# Website:          http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Provides:         fastqc
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/fastqc:v0.11.5 -f fastqc-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/fastqc:v0.11.5
# Run Cmd:          docker run --rm -ti biowardrobe2/fastqc:v0.11.5 fastqc --version
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing fastqc

ENV VERSION 0.11.5
ENV NAME fastqc
ENV URL "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/${NAME}_v${VERSION}.zip"

RUN wget -q -O /opt/${NAME}_v${VERSION}.zip $URL && \
    cd /opt && \
    unzip ./${NAME}_v${VERSION}.zip &&\
    chmod 777 /opt/FastQC/${NAME} && \
    ln -s /opt/FastQC/${NAME} /usr/local/bin/ &&\
    rm /opt/${NAME}_v${VERSION}.zip
