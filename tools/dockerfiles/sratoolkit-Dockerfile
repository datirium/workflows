#################################################################
# Dockerfile
#
# Software:         sratoolkit
# Software Version: 2.8.2-1
# Description:      sratoolkit image for BioWardrobe
# Website:          https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
# Provides:         sratoolkit
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/sratoolkit:v2.8.2-1 -f sratoolkit-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/sratoolkit:v2.8.2-1
# Run Cmd:          docker run --rm -ti biowardrobe2/sratoolkit:v2.8.2-1 fastq-dump --help
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing sratoolkit

ENV VERSION 2.8.2-1
ENV NAME sratoolkit
ENV URL "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
    cp -r sratoolkit.${VERSION}-ubuntu64/bin/* /usr/local/bin && \
    cp -r sratoolkit.${VERSION}-ubuntu64/schema /usr/local && \
    rm -rf sratoolkit.${VERSION}-ubuntu64
