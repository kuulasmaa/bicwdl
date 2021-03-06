# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="0.1"
LABEL description="UEF Bioinformatics Center / bioinformatics toolkit for GWAS analysis"
LABEL maintainer="Teemu Kuulasmaa <teemu.kuulasmaa@uef.fi>"

# Environment
ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/bin/fastqc/:/usr/local/bin/gatk/:${PATH}"
ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV HOME /root

# Define default command.
CMD ["bash"]

# System
RUN apt-get -q update && \
    apt-get -q -y install --no-install-recommends libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev \
    libncurses5 zlib1g libbz2-1.0 liblzma5 libssl1.1 libcurl4 \
    build-essential automake autoconf cmake git curl unzip dos2unix datamash ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# bcftools
RUN curl -SsOL https://github.com/samtools/bcftools/releases/download/1.10/bcftools-1.10.tar.bz2 && \
    tar jxf bcftools-1.10.tar.bz2
WORKDIR /bcftools-1.10
RUN ./configure && make && make install
WORKDIR /

# htslib
RUN curl -SsOL https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 && \
    tar jxf htslib-1.10.tar.bz2
WORKDIR /htslib-1.10
RUN ./configure && make && make install
WORKDIR /
RUN rm -rf bcftools-1.10* htslib-1.10*

# plink1.9, plink2
RUN curl -SsOL http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200121.zip && unzip plink_linux_x86_64_20200121.zip -d plink
RUN rm plink_linux_x86_64_20200121.zip
RUN cp plink/plink /usr/local/bin/
RUN rm -rf plink

RUN curl -SsOL http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20200208.zip && unzip plink2_linux_x86_64_20200208.zip -d plink2
RUN rm plink2_linux_x86_64_20200208.zip
RUN cp plink2/plink2 /usr/local/bin/
RUN rm -rf plink2

# clean-up
#  apt-utils apt-transport-https software-properties-common
RUN apt-get -q -y purge build-essential curl automake autoconf cmake libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev
RUN apt-get -q -y autoremove
RUN apt-get -q -y clean
RUN rm -rf /var/lib/apt/lists/* /var/lib/apt/extended_states \
  /var/cache/* /var/lib/dpkg/status /var/lib/dpkg/status-old \
  /etc/ld.so.cache /var/log/faillog /var/log/lastlog \
  /var/lib/apt/extended_states

# Define working directory.
WORKDIR /
