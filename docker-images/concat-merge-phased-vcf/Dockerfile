FROM debian:buster
  
LABEL \
  author="Dustin Miller" \
  description="miniconda2, bcftools" \
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
  # Install samtools
  && conda install -c bioconda bcftools=1.9 \
  && conda install -c bioconda samtools=1.9

ENV \
  PATH=${PATH}:/root/miniconda2/bin/