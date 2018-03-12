#################################################################
# Dockerfile
#
# Software:         trimmomatic
# Software Version: 0.35+dfsg-1
# Description:      trimmomatic image for BioWardrobe
# Website:          http://www.usadellab.org/cms/?page=trimmomatic
# Provides:         trimmomatic
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/trimmomatic:v0.35 -f trimmomatic-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/trimmomatic:v0.35
# Run Cmd:          docker run --rm biowardrobe2/trimmomatic:v0.35
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing trimmomatic

ENV VERSION 0.35+dfsg-1
ENV NAME trimmomatic

RUN apt-get update && apt-get install -y \
  ${NAME}="${VERSION}"