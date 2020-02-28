FROM debian:buster
  
LABEL \
  author="Dustin Miller" \
  description="Install miniconda2, GATK4, and reference files" \
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
  && conda config --append channels anaconda \
  # Download Reference Files
  && wget --no-check-certificate https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d9f54d2a7bc73000ee99fd6/?zip= -O /references.zip \
  && unzip /references.zip -d /references \
  && rm /references.zip \
  && gzip -d /references/*.gz \
  # Install GATK4
  && conda install -c bioconda gatk4=4.0.5.1 \
  # Install bcftools
  && conda install -c bioconda bcftools=1.9 \
  # Install samtools
  && conda install -c bioconda samtools=1.9
ENV \
  PATH=${PATH}:/root/miniconda2/bin/