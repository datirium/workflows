#################################################################
# Dockerfile
#
# Software:         Homer
# Software Version: v4.10
# Description:      Suite of tools for Motif Discovery and next-gen sequencing analysis
# Website:          http://homer.ucsd.edu/homer/
# Provides:         Homer basic installation
# Base Image:       biowardrobe2/deseq:v0.0.1
# Build Cmd:        docker build --rm -t biowardrobe2/homer:v0.0.2 -f homer-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/homer:v0.0.2
# Run Cmd:          docker run --rm -ti biowardrobe2/homer:v0.0.2
#################################################################
# v0.0.2
# Patched annotatePeaks.pl script to export all regions
#   remove line #3094 - `next if (!exists($ghistData{$peakID}));`
#
# v0.0.1
# The basic Homer installation does not contain any sequence data.
# To get a list of available packages: perl configureHomer.pl -list
#################################################################

### Base Image
FROM biowardrobe2/deseq:v0.0.1
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV VER_HOMER v4.10
ENV URL_HOMER "http://homer.ucsd.edu/homer/configureHomer.pl"
ENV DIR_HOMER "/opt/homer"

ENV PATH "$PATH:$DIR_HOMER/bin"

COPY ./scripts/install_from_bioconductor.R /tmp/install_from_bioconductor.R

### Installing edgeR package from Bioconductor
RUN Rscript /tmp/install_from_bioconductor.R edgeR && \

### Installing Homer
    mkdir ${DIR_HOMER} && \
    cd ${DIR_HOMER} && \
    wget -q ${URL_HOMER} && \
    perl configureHomer.pl -install homer -version ${VER_HOMER} && \

### Patch annotatePeaks.pl script
    cat /opt/homer/bin/annotatePeaks.pl | awk '{if (NR != 3094) print $0}' > /opt/homer/bin/annotatePeaks_updated.pl && \
    mv /opt/homer/bin/annotatePeaks_updated.pl /opt/homer/bin/annotatePeaks.pl && \
    chmod +x /opt/homer/bin/annotatePeaks.pl && \

### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
