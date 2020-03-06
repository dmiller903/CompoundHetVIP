FROM debian:buster
  
LABEL \
  author="Dustin Miller" \
  description="Tools needed for Compound Heterozygous Variant Identification" \
  maintainer="dustin_miller@byu.edu"

RUN \
# Install dependencies and other useful programs
apt-get update \
  && apt-get install -y --no-install-recommends wget=1.20.1-1.1 \
  build-essential=12.6 \
  pkg-config=0.29-6 \
  python3=3.7.3-1 \
  vim=2:8.1.0875-5 \
  unzip=6.0-23+deb10u1 \
  tar=1.30+dfsg-6 \
  zlib1g-dev=1:1.2.11.dfsg-1 \
  autotools-dev=20180224.1 \
  automake=1:1.16.1-4 \
  autoconf=2.69-11 \
  autogen=1:5.18.12-4 \
  # Install miniconda2 and config channels
  && wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda2-4.7.10-Linux-x86_64.sh \
  && bash Miniconda2-4.7.10-Linux-x86_64.sh -b -p \
  && export PATH="/root/miniconda2/bin:$PATH" \
  && /root/miniconda2/bin/conda init \
  && conda config --append channels bioconda \
  && conda config --append channels cyclus \
  && conda config --append channels conda-forge \
  # Install xz
  && conda install -c conda-forge xz=5.2.4 \
  # Install Java
  && conda install -c bioconda java-jdk=8.0.112 \
  # Install Plink2
  && wget --no-check-certificate https://github.com/dmiller903/ch-pipeline/raw/master/plink2_linux_x86_64_20190928.zip \
  && unzip /plink2_linux_x86_64_20190928.zip \
  # Install Picard
  && conda install -c bioconda picard=2.21.1 \
  # Install GATK4
  && conda install -c bioconda gatk4=4.0.5.1 \
  # Install bcftools
  && conda install -c bioconda bcftools=1.9 \
  # Install samtools
  && conda install -c bioconda samtools=1.9 \
  # Install shapeit
  && wget --no-check-certificate https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz \
  && tar -zxvf shapeit.v2.r904.glibcv2.17.linux.tar.gz \
  # Download Beagle Java file
  && wget --no-check-certificate http://faculty.washington.edu/browning/beagle/beagle.25Nov19.28d.jar \
  # Download and Extract Eagle executable
  && wget --no-check-certificate https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz \
  && tar xzf Eagle_v2.4.1.tar.gz \
  && rm Eagle_v2.4.1.tar.gz \
  && echo "eagle='./Eagle_v2.4.1/eagle'" > .bash_profile \
  # Install vt
  && conda install -c bioconda vt=2015.11.10 \
  # Install grabix (needed for GEMINI)
  && conda install -c bioconda grabix=0.1.8 \
  # Install tabix (needed for GEMINI)
  && conda install -c bioconda tabix=0.2.6 \
  # Install snpEff
  && wget --no-check-certificate https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip \
  && unzip /snpEff_v4_3t_core.zip \
  && mkdir -p /snpEff/./data/GRCh37.75 \
  # Install gemini without annotation files or cadd scores
  && wget --no-check-certificate https://raw.githubusercontent.com/dmiller903/gemini/master/gemini/scripts/gemini_install.py \
  && python gemini_install.py --nodata /usr/local /usr/local/share/gemini
ENV \
  PATH=${PATH}:/root/miniconda2/bin/:/root/miniconda2/share/snpeff-4.3.1t-2:/root/miniconda2/bin:/usr/local/bin/gemini