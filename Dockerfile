FROM kspacekelvin/fire-python-devcon

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git make curl sudo xvfb bzip2 gcc g++ python3-pip \
                    libglib2.0-0 libgl1 libxrender1 libxkbcommon-x11-0 libdbus-1-3 && \
    apt-get clean

# Clone the latest version of python-ismrmrd-server
RUN mkdir -p /opt/code/python-ismrmrd-server && \
    git clone https://github.com/kspaceKelvin/python-ismrmrd-server.git /opt/code/python-ismrmrd-server

# Install SCT
RUN mkdir -p /opt/code/spinal-cord-toolbox && \
    git clone --branch 7.0 --depth 1 https://github.com/spinalcordtoolbox/spinalcordtoolbox.git /opt/code/spinal-cord-toolbox && \
    cd /opt/code/spinal-cord-toolbox && \
    ./install_sct -y

# Install ST
RUN mkdir -p /opt/code/shimming-toolbox && \
    git clone --branch 1.2 --depth 1 https://github.com/shimming-toolbox/shimming-toolbox.git /opt/code/shimming-toolbox && \
    cd /opt/code/shimming-toolbox && \
    make install PLUGIN=false

# Install prelude
RUN mkdir /opt/code/prelude && \
    curl -o /opt/code/prelude/prelude -JL https://github.com/shimming-toolbox/binaries/raw/master/prelude && \
    sudo install /opt/code/prelude/prelude /usr/local/bin

# Activate ST environment
SHELL ["/bin/bash", "-c"]
RUN echo 'source /opt/code/shimming-toolbox/python/etc/profile.d/conda.sh' >> ~/.bashrc && \
    echo 'conda activate /opt/code/shimming-toolbox/python' >> ~/.bashrc

# Set the starting directory so that code can use relative paths
WORKDIR /opt/code/shimming-toolbox

# Clean up to reduce image size
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    rm -rf /opt/code/spinal-cord-toolbox/.git /opt/code/shimming-toolbox/.git /opt/code/python-ismrmrd-server/.git && \
    rm -rf /root/.cache/pip && \
    rm -rf /var/cache/apt/archives && \
    find /opt/code -name "*.pyc" -delete && \
    find /opt/code -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true && \
    rm -rf /opt/code/*/docs /opt/code/*/examples /opt/code/*/tests && \
    rm -rf /opt/code/spinal-cord-toolbox/data/PAM50 /opt/code/spinal-cord-toolbox/data/template && \
    apt-get autoremove -y && \
    apt-get autoclean
