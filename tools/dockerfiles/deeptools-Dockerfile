#################################################################
# Dockerfile
#
# Software:         deepTools
# Software Version: v0.0.1
# Description:      Dockerized version of deepTools
# Website:          https://github.com/Barski-lab/workflows
# Provides:         deepTools 3.3.0
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/deeptools:v0.0.1 -f deeptools-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/deeptools:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/deeptools:v0.0.1
#################################################################


### Base Image
FROM biowardrobe2/scidap:v0.0.3
LABEL maintainer="misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


ENV VER_DEEPTOOLS 3.3.0

################## BEGIN INSTALLATION ######################
WORKDIR /tmp

### Installing deepTools
RUN pip install -U pip && \
    hash pip && \
    pip install setuptools && \
    pip install deeptools==${VER_DEEPTOOLS} && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
