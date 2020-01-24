FROM debian:buster
  
LABEL \
  author="Dustin Miller" \
  description="shapeit, HapCUT2, Beagle, and Eagle" \
  maintainer="dustin_miller@byu.edu"

RUN \
# Install dependencies and other useful programs
apt-get update \
  && apt-get install -y --no-install-recommends wget=1.20.1-1.1 \
  python3=3.7.3-1 \
  build-essential=12.6 \
  pkg-config=0.29-6 \
  vim=2:8.1.0875-5 \
  unzip \
  tar \
  zlib1g-dev=1:1.2.11.dfsg-1 \
  autotools-dev=20180224.1 \
  automake=1:1.16.1-4 \
  autoconf=2.69-11 \
  autogen=1:5.18.12-4 \
  # Install miniconda2
  && wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda2-4.7.10-Linux-x86_64.sh \
  && bash Miniconda2-4.7.10-Linux-x86_64.sh -b -p \
  && export PATH="/root/miniconda2/bin:$PATH" \
  && /root/miniconda2/bin/conda init \
  && conda config --append channels bioconda \
  && conda config --append channels conda-forge \
  && conda config --append channels anaconda \
  && conda config --append channels cyclus \
  # Install other dependencies through conda
  && conda install -c conda-forge xz \
  && conda install -c cyclus java-jdk \
  # Install git
  && conda install -c anaconda git=2.20.1 \
  # Download Reference Files
  && wget --no-check-certificate https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz -P /references \
  && cd /references \
  && tar zxvf /references/1000GP_Phase3.tgz \
  && cd / \
  # Install shapeit
  && wget --no-check-certificate https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz \
  && tar -zxvf shapeit.v2.r904.glibcv2.17.linux.tar.gz \
  # Install HapCUT2
  && git clone https://github.com/dmiller903/HapCUT2.git \
  && cd /HapCUT2 \
  && make \
  && make install \
  && cd / \
  # Download Beagle Java file
  && wget --no-check-certificate http://faculty.washington.edu/browning/beagle/beagle.25Nov19.28d.jar \
  # Download and Extract Eagle executable
  && wget --no-check-certificate https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz \
  && tar xzf Eagle_v2.4.1.tar.gz \
  && rm Eagle_v2.4.1.tar.gz \
  && echo "eagle='./Eagle_v2.4.1/eagle'" > .bash_profile
ENV \
  PATH=${PATH}:/root/miniconda2/bin
