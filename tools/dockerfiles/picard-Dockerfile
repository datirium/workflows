#################################################################
# Dockerfile
#
# Software:         picard
# Software Version: 2.9.0
# Description:      picard image for BioWardrobe
# Website:          https://broadinstitute.github.io/picard/
# Provides:         picard
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/picard:v2.9.0 -f picard-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/picard:v2.9.0
# Run Cmd:          docker run --rm -ti biowardrobe2/picard:v2.9.0 java -jar /usr/local/bin/picard.jar
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing picard

ENV VERSION 2.9.0
ENV NAME picard
ENV URL "https://github.com/broadinstitute/picard/releases/download/${VERSION}/picard.jar"

RUN wget -q $URL && \
    mv picard.jar /usr/local/bin/
