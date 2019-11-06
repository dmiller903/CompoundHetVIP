FROM debian:buster

LABEL \
  author="Dustin Miller" \
  description="GEMINI, snpEff, and dependencies" \
  maintainer="dustin_miller@byu.edu"

RUN \
# Install dependencies and other useful programs
apt-get update \
  && apt-get install -y --no-install-recommends wget=1.20.1-1.1 \
  python3=3.7.3-1 \
  #build-essential=12.6 \
  #pkg-config=0.29-6 \
  vim=2:8.1.0875-5 \
  unzip=6.0-23+deb10u1 \
  #zlib1g-dev=1:1.2.11.dfsg-1 \
  #autotools-dev=20180224.1 \
  #automake=1:1.16.1-4 \
  #autoconf=2.69-11 \
  #autogen=1:5.18.12-4 \
  # Install miniconda2
  && wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda2-4.7.10-Linux-x86_64.sh \
  && bash Miniconda2-4.7.10-Linux-x86_64.sh -b -p \
  && export PATH="/root/miniconda2/bin:$PATH" \
  && /root/miniconda2/bin/conda init \
  && conda config --append channels bioconda \
  && conda config --append channels conda-forge \
  && conda config --append channels anaconda \
  # Install other dependencies through conda
  && conda install -c conda-forge xz \
  && conda install -c bioconda grabix \
  && conda install -c bioconda tabix \
  # Install snpEff
  && conda install -c bioconda snpeff=4.3.1t \
  # Download Reference Files
  && wget --no-check-certificate https://phoenixnap.dl.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh37.75.zip -P /tmp \
  && unzip /tmp/snpEff_v4_3_GRCh37.75.zip -d /root/miniconda2/share/snpeff-4.3.1t-2/ \
  # Install gemini
  && wget --no-check-certificate https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py \
  && python gemini_install.py /usr/local /usr/local/share/gemini \
  && export PATH=$PATH:/usr/local/gemini/bin

ENV \
  PATH=${PATH}:/root/miniconda2/bin:/root/miniconda2/share/snpeff-4.3.1t-2