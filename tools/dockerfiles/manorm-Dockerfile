#################################################################
# Dockerfile
#
# Software:         BioWardrobe MAnorm script
# Software Version: v0.0.2
# Description:      Dockerized version of MAnorm Python Package
# Website:          https://github.com/Barski-lab/workflows
# Provides:         MAnorm v1.2.0
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/manorm:v0.0.2 -f manorm-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/manorm:v0.0.2
# Run Cmd:          docker run --rm -ti biowardrobe2/manorm:v0.0.2
#################################################################
#
# v0.0.2
# MAnorm Python Package v1.2.0 from the latest available commit (Pypi version is outdated v1.1.4)
# https://github.com/shao-lab/MAnorm/tree/82f3db7909caa3b58da9ba761524f32fb0d6b071
#
# v0.0.1
# MAnorm.R from BioWardrobe
#
#################################################################


### Base Image
FROM biowardrobe2/scidap:v0.0.3
LABEL maintainer="misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


ENV VER_MANORM 82f3db7909caa3b58da9ba761524f32fb0d6b071
ENV URL_MANORM git+git://github.com/shao-lab/MAnorm.git@${VER_MANORM}

################## BEGIN INSTALLATION ######################
WORKDIR /tmp

### Installing MAnorm
RUN pip install -U pip && \
    hash pip && \
    pip install setuptools && \
    pip install ${URL_MANORM} && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
