#################################################################
# Dockerfile
# Software:         bismark
# Software Version: v0.0.2
# Description:      bismark image for BioWardrobe2
# Website:          https://www.bioinformatics.babraham.ac.uk/projects/bismark/
# Provides:         bismark 0.22.1
#                   bowtie 1.2.0
#                   bowtie2 2.3.5.1
#                   htslib 1.9
#                   samtools 1.9
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/bismark:v0.0.2 -f bismark-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/bismark:v0.0.2
# Run Cmd:          docker run --rm -ti biowardrobe2/bismark:v0.0.2
#################################################################
# 
# v0.0.2
#   Updates:  bismark 0.19.0  --> 0.22.1
#             bowtie2 2.3.0   --> 2.3.5.1
#             htslib 1.4      --> 1.9
#             samtools 1.4    --> 1.9
#
#################################################################




### Base Image
FROM biowardrobe2/scidap:v0.0.3
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV VERSION_BISMARK 0.22.1
ENV URL_BISMARK "https://github.com/FelixKrueger/Bismark/archive/${VERSION_BISMARK}.tar.gz"

ENV VERSION_BOWTIE 1.2.0
ENV URL_BOWTIE "https://github.com/BenLangmead/bowtie/archive/v${VERSION_BOWTIE}.tar.gz"

ENV VERSION_BOWTIE2 2.3.5.1
ENV URL_BOWTIE2 "https://github.com/BenLangmead/bowtie2/archive/v${VERSION_BOWTIE2}.tar.gz"

ENV VERSION_HTSLIB 1.9
ENV URL_HTSLIB "https://github.com/samtools/htslib/releases/download/${VERSION_HTSLIB}/htslib-${VERSION_HTSLIB}.tar.bz2"

ENV VERSION_SAMTOOLS 1.9
ENV URL_SAMTOOLS "https://github.com/samtools/samtools/releases/download/${VERSION_SAMTOOLS}/samtools-${VERSION_SAMTOOLS}.tar.bz2"

### Installing dependencies
RUN apt-get update && \
    apt-get install -y libtbb-dev libgd-dev && \
    apt-get install -y libncurses5-dev libbz2-dev liblzma-dev && \
### Installing bismark
    wget -q -O /opt/Bismark-${VERSION_BISMARK}.tar.gz ${URL_BISMARK} && \
    cd /opt && \
    tar xzf ./Bismark-${VERSION_BISMARK}.tar.gz && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/NOMe_filtering /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/bam2nuc /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark2bedGraph /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark2report /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark2summary /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark_genome_preparation /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark_methylation_extractor /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/coverage2cytosine /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/deduplicate_bismark /usr/local/bin/ && \
    ln -s /opt/Bismark-${VERSION_BISMARK}/filter_non_conversion /usr/local/bin/ && \
    rm /opt/Bismark-${VERSION_BISMARK}.tar.gz && \
### Installing bowtie
    wget -q -O - $URL_BOWTIE | tar -zxv && \
    cd bowtie-${VERSION_BOWTIE} && \
    make -j 4 && \
    make install && \
    cd .. && \
    rm -rf ./bowtie-${VERSION_BOWTIE} && \
### Installing bowtie2
    wget -q -O - $URL_BOWTIE2 | tar -zxv && \
    cd bowtie2-${VERSION_BOWTIE2} && \
    make -j 4 && \
    make install && \
    cd .. && \
    rm -rf ./bowtie2-${VERSION_BOWTIE2} && \
### Installing samtools
    wget -q $URL_HTSLIB && \
    bzip2 -d htslib-${VERSION_HTSLIB}.tar.bz2 && \
    tar -xf htslib-${VERSION_HTSLIB}.tar && \
    cd htslib-${VERSION_HTSLIB} && \
    ./configure && \
    make -j 4 && \
    make install && \
    cd .. && \
    rm -rf ./htslib-${VERSION_HTSLIB} && \
    rm -rf ./htslib-${VERSION_HTSLIB}.tar && \
    wget -q $URL_SAMTOOLS && \
    bzip2 -d samtools-${VERSION_SAMTOOLS}.tar.bz2 && \
    tar -xf samtools-${VERSION_SAMTOOLS}.tar && \
    cd samtools-${VERSION_SAMTOOLS} && \
    ./configure && \
    make -j 4 && \
    make install &&\
    cd .. && \
    rm -rf ./samtools-${VERSION_SAMTOOLS} && \
    rm -rf ./samtools-${VERSION_SAMTOOLS}.tar && \
### Installing required Perl modules
    perl -MCPAN -e 'my $c = "CPAN::HandleConfig"; $c->load(doit => 1, autoconfig => 1); $c->edit(prerequisites_policy => "follow"); $c->edit(build_requires_install_policy => "yes"); $c->commit' && \
    cpan App::cpanminus && \
    cpanm GD::Graph && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true