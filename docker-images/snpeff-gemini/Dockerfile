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
  vim=2:8.1.0875-5 \
  unzip=6.0-23+deb10u1 \
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
  # Install gemini with cadd scores
  && wget --no-check-certificate https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py \
  && python gemini_install.py /usr/local /usr/local/share/gemini \
  && export PATH=$PATH:/usr/local/gemini/bin \
  && gemini update --dataonly --extra cadd_score

ENV \
  PATH=${PATH}:/root/miniconda2/bin:/root/miniconda2/share/snpeff-4.3.1t-2