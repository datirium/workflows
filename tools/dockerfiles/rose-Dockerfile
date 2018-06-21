#################################################################
# Dockerfile
#
# Software:         ROSE
# Software Version: v0.0.2
# Description:      ROSE image for BioWardrobe
# Website:          https://github.com/Barski-lab/ROSE/
# Provides:         ROSE (RANK ORDERING OF SUPER-ENHANCERS)
# Base Image:       biowardrobe2/scidap-addons:v0.0.6
# Build Cmd:        docker build --rm -t biowardrobe2/rose:v0.0.2 -f rose-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/rose:v0.0.2
# Run Cmd:          docker run --rm -ti biowardrobe2/rose:v0.0.2 /bin/bash
#################################################################


### Base Image
FROM biowardrobe2/scidap-addons:v0.0.6
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV VERSION_ROSE 0.0.2
ENV URL_ROSE "https://github.com/Barski-lab/ROSE/archive/v${VERSION_ROSE}.tar.gz"

COPY ./scripts/rose.py /usr/local/bin/ROSE_main
COPY ./scripts/makegff.R /usr/local/bin/makegff

### Installing ROSE
RUN mkdir /opt/rose && \
    cd /opt/rose && \
    wget -q -O - $URL_ROSE | tar -zxv --strip-components=1 &&\
    chmod +x /usr/local/bin/ROSE_main && \

### Installing makegff R script
    chmod +x /usr/local/bin/makegff && \

### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
