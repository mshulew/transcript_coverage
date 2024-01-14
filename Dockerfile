FROM ubuntu:18.04
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    perl \
    default-jdk \
    wget \
    samtools \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    pigz \
    python \
    python3.6 \
    python3-pip
    
######### Python 3.9 #############   
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        software-properties-common \
    && add-apt-repository -y ppa:deadsnakes \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3.9

######### Subread Setup #############
ENV SUBREAD_VERSION 1.6.4 
RUN curl -SLO https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz \
    && tar -zxvf subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz --directory /opt/ \
    && rm subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz
ENV PATH /opt/subread-${SUBREAD_VERSION}-Linux-x86_64/bin/:$PATH
######### End Subread Setup ##########

######### bedtools Setup #############
ENV BEDTOOLS_VERSION 2.30.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
&& tar -zxvf bedtools-${BEDTOOLS_VERSION}.tar.gz && rm bedtools-${BEDTOOLS_VERSION}.tar.gz
RUN make -C bedtools2
RUN mv bedtools2/bin /opt/bedtools
RUN rm -r bedtools2
ENV PATH /opt/bedtools/:$PATH
######### End bedtools Setup ##########

######### Sambamba Setup #############
ENV SAMBAMBA_VERSION 0.8.0
RUN mkdir -p /opt/sambamba
RUN curl -SLO https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba-${SAMBAMBA_VERSION}-linux-amd64-static.gz \
    && unpigz sambamba-${SAMBAMBA_VERSION}-linux-amd64-static.gz && mv sambamba-${SAMBAMBA_VERSION}-linux-amd64-static /opt/sambamba/sambamba
RUN chmod 775 /opt/sambamba/sambamba
ENV PATH /opt/sambamba/:$PATH
######### End Sambamba Setup #########

WORKDIR /opt/biorad

COPY . .
