# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="0.1"
LABEL description="UEF Bioinformatics Center / bioinformatics toolkit for NGS analysis"
LABEL maintainer="Teemu Kuulasmaa <teemu.kuulasmaa@uef.fi>"

# Environment
ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/bin/fastqc/:/usr/local/bin/gatk/:${PATH}"
ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV HOME /root
ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/

# Define default command.
CMD ["bash"]

# System
RUN apt-get -q update && \
    apt-get -q -y install --no-install-recommends libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev \
    libncurses5 zlib1g libbz2-1.0 liblzma5 libssl1.1 libcurl4 openjdk-11-jre-headless \
    build-essential automake autoconf cmake pkg-config libgtextutils0v5 libgtextutils-dev file git curl unzip ca-certificates dos2unix datamash python3-setuptools python3-wheel python3-pip && \
    ln -sf /usr/bin/python3 /usr/bin/python && rm -rf /var/lib/apt/lists/*

# fastx-toolkit 
# > compile error: https://github.com/agordon/fastx_toolkit/issues/21
#RUN curl -SsOL https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2 && \
#    tar jxf fastx_toolkit-0.0.14.tar.bz2
#WORKDIR /fastx_toolkit-0.0.14
#RUN ./configure && make && make install
#WORKDIR /
#RUN rm -rf fastx_toolkit-0.0.14.tar.bz2
RUN apt-get install fastx-toolkit
  
# bwa
RUN curl -SsOL https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar jxf bwa-0.7.17.tar.bz2
WORKDIR /bwa-0.7.17
RUN make && cp bwa /usr/local/bin/
WORKDIR /
RUN rm -rf bwa-0.7.17*

# samtools
RUN curl -SsOL https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar jxf samtools-1.10.tar.bz2
WORKDIR /samtools-1.10
RUN ./configure && make && make install
WORKDIR /
RUN rm -rf samtools-1.10*

# bcftools
RUN curl -SsOL https://github.com/samtools/bcftools/releases/download/1.10/bcftools-1.10.tar.bz2 && \
    tar jxf bcftools-1.10.tar.bz2
WORKDIR /bcftools-1.10
RUN ./configure && make && make install
WORKDIR /
RUN rm -rf bcftools-1.10*

# htslib
RUN curl -SsOL https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 && \
    tar jxf htslib-1.10.tar.bz2
WORKDIR /htslib-1.10
RUN ./configure && make && make install
WORKDIR /
RUN rm -rf htslib-1.10*

# kallisto
RUN curl -SsOL https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar zxf kallisto_linux-v0.46.1.tar.gz && \
    cp kallisto/kallisto /usr/local/bin/ && rm -rf kallisto*

# fastqc
RUN curl -SsOL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    rm fastqc_v0.11.9.zip && rm -rf FastQC/Help && \
    mv FastQC /usr/local/bin/ && chmod 755 /usr/local/bin/FastQC/fastqc && \
    ln -sf /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc

# gatk
#RUN curl -SsOL https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip && \
#    unzip gatk-4.1.4.1.zip && \
#    rm gatk-4.1.4.1.zip && rm -rf gatk-4.1.4.1/gatkdoc && \
#    rm -rf gatk-4.1.4.1/scripts && \
#    mv gatk-4.1.4.1 /usr/local/bin/GATK && chmod 755 /usr/local/bin/GATK/gatk && \
#    ln -sf /usr/local/bin/GATK/gatk /usr/local/bin/gatk

# bowtie2
RUN curl -SsOL https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.5.1-linux-x86_64.zip -d bowtie2 && \
    rm bowtie2-2.3.5.1-linux-x86_64.zip && \
    mv bowtie2 /usr/local/bin/ && \
    ln -sf /usr/local/bin/bowtie2/bowtie2 /usr/local/bin/bowtie2

# Trimmomatic
RUN curl -SsOL ttp://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    rm Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /usr/local/bin/ && \
    ln -sf /usr/local/bin/Trimmomatic-0.39/Trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar

# cutadapt, multiqc
RUN pip3 install cutadapt==2.8 multiqc==1.8

# clean-up
RUN apt-get -q -y purge build-essential curl automake autoconf cmake libgtextutils-dev libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev 
RUN apt-get -q -y autoremove
RUN apt-get -q -y clean
RUN rm -rf /var/lib/apt/lists/* /var/lib/apt/extended_states \
  /var/cache/* /var/lib/dpkg/status /var/lib/dpkg/status-old \
  /etc/ld.so.cache /var/log/faillog /var/log/lastlog \
  /var/lib/apt/extended_states

# Define working directory.
WORKDIR /



