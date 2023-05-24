FROM ghcr.io/iqbal-lab-org/clockwork:v0.11.3

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN apt-get update && apt-get install -y zlib1g-dev wget libbz2-dev liblzma-dev bzip2 make gcc libncurses5-dev \
					 automake build-essential cmake curl gawk git gnuplot graphviz libarchive-dev libcurl4-gnutls-dev \
					 libhts-dev libncursesw5-dev libvcflib-tools libssl-dev pkg-config python3-dev python3-pip \
					 python3-setuptools r-base-core rsync unzip tabix libssl-dev openjdk-8-jdk

get -q https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip
RUN printf '#!/bin/bash \
\njava -jar /Trimmomatic-0.39/trimmomatic-0.39.jar "$@"\n' > /usr/local/bin/trimmomatic && chmod u+x /usr/local/bin/trimmomatic

RUN wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && unzip snpEff_latest_core.zip
RUN printf '#!/bin/bash \
\njava -jar /snpEff/snpEff.jar $@"\n' > /usr/local/bin/snpEff && chmod u+x /usr/local/bin/snpEff

RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.2.4.1/gatk-4.2.4.1.zip && unzip gatk-4.2.4.1.zip

RUN apt-get install -y bwa bedtools 

ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64
#RUN curl -fsSL get.nextflow.io | bash
#RUN mv nextflow /usr/local/bin

#RUN cd /varpipe_wgs/tools && wget -qO- https://get.nextflow.io | bash && chmod +x nextflow
#RUN cp -R /clockwork/nextflow /varpipe_wgs/tools/clockwork-0.11.3

COPY ./tools /varpipe_wgs/tools
WORKDIR /varpipe_wgs/data

ENV PATH=$PATH:/varpipe_wgs/tools:/gatk-4.2.4.1

CMD /bin/bash
