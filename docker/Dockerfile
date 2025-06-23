# Preamble
FROM ubuntu:latest
#FROM --platform=linux/amd64 ubuntu:latest

WORKDIR /
SHELL ["/bin/bash", "-c"]

# Install dependencies
RUN apt-get update && \
  apt-get install --no-install-recommends -y \
  curl git g++ ca-certificates && \

# Install Miniforge
  curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
  bash Miniforge3-$(uname)-$(uname -m).sh -b -p /opt/miniforge3/ && \

# Install Exorcise
  git clone https://github.com/SimonLammmm/exorcise && \
  /opt/miniforge3/bin/mamba env create -f /exorcise/env/exorcise.yaml -y python=3.12 && \
  export PATH=$PATH:/exorcise/bin && \
# Patch the Exorcise main script to use the blat installation in Miniforge
  sed -i 's/blat_command <- "blat"/blat_command <- "\/opt\/miniforge3\/envs\/exorcise\/bin\/blat"/g' /exorcise/bin/exorcise.R && \
  sed -i 's/twoBitToFa_command <- "twoBitToFa"/twoBitToFa_command <- "\/opt\/miniforge3\/envs\/exorcise\/bin\/twoBitToFa"/g' /exorcise/bin/exorcise.R && \

# Install crispr_tools
  git clone https://github.com/SimonLammmm/crispr_tools.git && \
# Patch the crispr_pipeline script to use the MAGeCK installation in Miniforge
  sed -i 's/mageck test/\/opt\/miniforge3\/envs\/crispr_tools\/bin\/mageck test/g' /crispr_tools/crispr_tools/crispr_pipeline.py && \
  sed -i 's/which mageck/which \/opt\/miniforge3\/envs\/crispr_tools\/bin\/mageck/g' /crispr_tools/crispr_tools/crispr_pipeline.py && \
  sed -i 's/"mageck", "-v"/"\/opt\/miniforge3\/envs\/crispr_tools\/bin\/mageck", "-v"/g' /crispr_tools/crispr_tools/crispr_pipeline.py && \
  /opt/miniforge3/bin/mamba env create -f /crispr_tools/env/crispr_tools.yaml -y && \
  /opt/miniforge3/envs/crispr_tools/bin/pip install /crispr_tools/. && \
# Patch the MAGeCK script to use RRA installation in Miniforge
  sed -i "s/rrapath='RRA'/rrapath='\/opt\/miniforge3\/envs\/crispr_tools\/bin\/RRA'/g" /opt/miniforge3/envs/crispr_tools/lib/python3.12/site-packages/mageck/crisprFunction.py && \

# Install crispr_screen_viewer
  git clone https://github.com/johncthomas/crispr_screen_viewer.git && \
  /opt/miniforge3/bin/mamba create -n crispr_screen_viewer python -y && \
  /opt/miniforge3/envs/crispr_screen_viewer/bin/pip install /crispr_screen_viewer/. scipy==1.15.3 && \

# Clean up
  /opt/miniforge3/bin/mamba clean --all -y && \
  rm -rf /opt/miniforge3/pkgs && \
  rm /Miniforge*.sh && \
  rm -rf /exorcise/data && \
  rm -rf /exorcise/env && \
  rm -rf /exorcise/example && \
  rm -rf /exorcise/manuscript && \
  rm -rf /crispr_tools/det && \
  rm -rf /crispr_tools/doc && \
  rm -rf /crispr_tools/env && \
  rm -rf /crispr_tools/files && \

# Install entrypoint
  chmod +x /exorcise/docker/entrypoint.sh
ENTRYPOINT ["/exorcise/docker/entrypoint.sh"]

# Mount working directory
VOLUME ["/data/"]
WORKDIR /data/
