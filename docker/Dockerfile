FROM rocker/r-ver:latest

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    git \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error \
    -n 4 \
    data.table \
    R.utils \
    magrittr \
    gridExtra \
    devtools \
    stringr \
    remotes
 
RUN R -e "remotes::install_github('twillis209/pidPipelineCode')"