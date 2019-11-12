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
  && conda install -c conda-forge xz=5.2.4 \
  && conda install -c bioconda grabix=0.1.8 \
  && conda install -c bioconda tabix=0.2.6 \
  && conda install -c bioconda java-jdk=8.0.112 \
  # Install snpEff and Reference Files
  && wget --no-check-certificate https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip \
  && unzip /snpEff_v4_3t_core.zip \
  && java -jar /snpEff/snpEff.jar download -v GRCh37.75 \
  && wget --no-check-certificate https://doc-08-cc-docs.googleusercontent.com/docs/securesc/ha0ro937gcuc7l7deffksulhg5h7mbp1/0ng6pgff7h5sdictk2hvvr0n84d18ogc/1573164000000/08741013876167739028/*/0B7Ms5xMSFMYlOTV5RllpRjNHU2s?e=download -O /snpEff/./data/dbNSFP2.9.txt.gz.tbi \
  && wget --no-check-certificate https://doc-0s-cc-docs.googleusercontent.com/docs/securesc/ha0ro937gcuc7l7deffksulhg5h7mbp1/ppsbqsfaf5fkfo9e7dvpkg2qno2iengt/1573164000000/08741013876167739028/*/0B7Ms5xMSFMYlSTY5dDJjcHVRZ3M?e=download -O /snpEff/./data/dbNSFP2.9.txt.gz \
  # Install gemini
  && wget --no-check-certificate https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py \
  && python gemini_install.py /usr/local /usr/local/share/gemini \
  && export PATH=$PATH:/usr/local/gemini/bin

ENV \
  PATH=${PATH}:/root/miniconda2/bin:/root/miniconda2/share/snpeff-4.3.1t-2