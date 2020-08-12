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
  git=1:2.20.1-2+deb10u3 \
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
  && wget --no-check-certificate "https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d950c8b4ab86b001af8b807?action=download&direct&version=1" -O /plink2.zip \
  && unzip /plink2.zip \
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
  # Install snpEff and create directory where snpEff files will eventually be stored
  && wget --no-check-certificate https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip \
  && unzip /snpEff_v4_3t_core.zip \
  && mkdir -p /snpEff/./data/GRCh37.75 \
  && mkdir -p /snpEff/./data/GRCh38.86 \
  # Install gemini without annotation files or cadd scores
  && wget --no-check-certificate https://raw.githubusercontent.com/dmiller903/gemini/master/gemini/scripts/gemini_install.py \
  && python gemini_install.py --nodata /usr/local /usr/local/share/gemini \
  # Clone CompoundHetVIP from github and place python scripts at root
  && git config --global http.sslverify "false" \
  && git clone https://github.com/dmiller903/CompoundHetVIP.git \
  && cp CompoundHetVIP/scripts/*.py / \
  && rm -r CompoundHetVIP \
  # Download GENCODE file
  && wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz -P / \
  # Download Gene Damage Index (GDI) Files
  && wget --no-check-certificate http://lab.rockefeller.edu/casanova/assets/file/add_GDI.py.txt -O /add_GDI_raw.py \
  && wget --no-check-certificate http://lab.rockefeller.edu/casanova/assets/file/GDI_full_10282015.txt -O /GDI_full.txt \
  # Install vcf2db
  && git clone https://github.com/dmiller903/vcf2db.git \
  && cd vcf2db \
  && conda install -y -c conda-forge snappy \
  && conda install -y -c conda-forge python-snappy \
  && conda install -y -c bioconda cyvcf2 \
  && conda install -y -c bioconda peddy \
  && pip install -r requirements.txt
ENV \
  PATH=${PATH}:/root/miniconda2/bin/:/root/miniconda2/share/snpeff-4.3.1t-2:/root/miniconda2/bin:/usr/local/bin/gemini