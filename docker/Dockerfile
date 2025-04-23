# Preamble
FROM ubuntu:latest
WORKDIR /
SHELL ["/bin/bash", "-c"]

# Install dependencies
RUN apt-get update && \
  apt-get install --no-install-recommends -y \
  curl git g++ ca-certificates

# Install Miniforge
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3-$(uname)-$(uname -m).sh -b

# Install Exorcise
RUN git clone https://github.com/SimonLammmm/exorcise
RUN /root/miniforge3/bin/mamba env create -f /exorcise/env/exorcise.yaml -y
RUN export PATH=$PATH:/exorcise/bin
# Patch the Exorcise main script to use the blat installation in Miniforge
RUN sed -i 's/blat_command <- "blat"/blat_command <- "\/root\/miniforge3\/envs\/exorcise\/bin\/blat"/g' /exorcise/bin/exorcise.R
RUN sed -i 's/twoBitToFa_command <- "twoBitToFa"/twoBitToFa_command <- "\/root\/miniforge3\/envs\/exorcise\/bin\/twoBitToFa"/g' /exorcise/bin/exorcise.R

# Install crispr_tools
RUN git clone https://github.com/SimonLammmm/crispr_tools.git
# Patch the crispr_pipeline script to use the MAGeCK installation in Miniforge
RUN sed -i 's/mageck test/\/root\/miniforge3\/envs\/crispr_tools\/bin\/mageck test/g' /crispr_tools/crispr_tools/crispr_pipeline.py
RUN sed -i 's/which mageck/which \/root\/miniforge3\/envs\/crispr_tools\/bin\/mageck/g' /crispr_tools/crispr_tools/crispr_pipeline.py
RUN sed -i 's/"mageck", "-v"/"\/root\/miniforge3\/envs\/crispr_tools\/bin\/mageck", "-v"/g' /crispr_tools/crispr_tools/crispr_pipeline.py
RUN /root/miniforge3/bin/mamba env create -f /crispr_tools/env/crispr_tools.yaml -y
RUN /root/miniforge3/envs/crispr_tools/bin/pip install /crispr_tools/.

# Install crispr_screen_viewer
RUN git clone https://github.com/johncthomas/crispr_screen_viewer.git
RUN /root/miniforge3/bin/mamba create -n crispr_screen_viewer python -y
RUN /root/miniforge3/envs/crispr_screen_viewer/bin/pip install /crispr_screen_viewer/.

# Clean up
RUN /root/miniforge3/bin/mamba clean --all -y

# Install entrypoint
RUN chmod +x /exorcise/docker/entrypoint.sh
ENTRYPOINT ["/exorcise/docker/entrypoint.sh"]

# Mount working directory
VOLUME ["/data/"]
WORKDIR /data/
