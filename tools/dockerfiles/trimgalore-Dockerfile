#################################################################
# Dockerfile
#
# Software:         TrimGalore [FastQC, Cutadapt]
# Software Version: 0.4.4      [0.11.5, 1.13]
# Description:      TrimGalore image for BioWardrobe
# Website:          https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
# Provides:         trim_galore, fastqc, cutadapt
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/trimgalore:v0.4.4 -f trimgalore-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/trimgalore:v0.4.4
# Run Cmd:          docker run --rm -ti biowardrobe2/trimgalore:v0.4.4
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing TrimGalore

ENV VERSION_FQC 0.11.5
ENV NAME_FQC fastqc
ENV URL_FQC "http://www.bioinformatics.babraham.ac.uk/projects/${NAME_FQC}/${NAME_FQC}_v${VERSION_FQC}.zip"

ENV VERSION_CADP 1.13
ENV NAME_CADP cutadapt

ENV VERSION_TRG 0.4.4
ENV NAME_TRG trim_galore
ENV URL_TRG "https://www.bioinformatics.babraham.ac.uk/projects/${NAME_TRG}/${NAME_TRG}_v${VERSION_TRG}.zip"


RUN wget -q -O /opt/${NAME_FQC}_v${VERSION_FQC}.zip ${URL_FQC} && \
    cd /opt && \
    unzip ./${NAME_FQC}_v${VERSION_FQC}.zip &&\
    chmod 777 /opt/FastQC/${NAME_FQC} && \
    ln -s /opt/FastQC/${NAME_FQC} /usr/local/bin/ &&\
    rm /opt/${NAME_FQC}_v${VERSION_FQC}.zip &&\
    pip install setuptools &&\
    pip install --upgrade pip &&\
    pip install ${NAME_CADP}==${VERSION_CADP} &&\
    mkdir /opt/${NAME_TRG} &&\
    wget --no-check-certificate -q -O /opt/${NAME_TRG}/${NAME_TRG}_v${VERSION_TRG}.zip ${URL_TRG} && \
    cd /opt/${NAME_TRG} && \
    unzip ./${NAME_TRG}_v${VERSION_TRG}.zip &&\
    chmod 777 /opt/${NAME_TRG}/${NAME_TRG} && \
    ln -s /opt/${NAME_TRG}/${NAME_TRG} /usr/local/bin/



