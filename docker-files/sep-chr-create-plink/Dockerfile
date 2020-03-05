FROM debian:buster
  
LABEL \
  author="Dustin Miller" \
  description="Install miniconda2, java, and plink2" \
  maintainer="dustin_miller@byu.edu"

RUN \
# Install dependencies and other useful programs
apt-get update \
  && apt-get install -y --no-install-recommends wget=1.20.1-1.1 \
  build-essential=12.6 \
  python3=3.7.3-1 \
  vim=2:8.1.0875-5 \
  unzip=6.0-23+deb10u1 \
  # Install miniconda2 and append channels
  && wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda2-4.7.10-Linux-x86_64.sh \
  && bash Miniconda2-4.7.10-Linux-x86_64.sh -b -p \
  && export PATH="/root/miniconda2/bin:$PATH" \
  && /root/miniconda2/bin/conda init \
  && conda config --append channels bioconda \
  && conda config --append channels anaconda \
  && conda config --append channels cyclus \
  # Install Java
  && conda install -c cyclus java-jdk \
  # Install Plink2
  && wget --no-check-certificate https://github.com/dmiller903/ch-pipeline/raw/master/plink2_linux_x86_64_20190928.zip \
  && unzip /plink2_linux_x86_64_20190928.zip

ENV \
  PATH=${PATH}:/root/miniconda2/bin/