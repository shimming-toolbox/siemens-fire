FROM kspacekelvin/fire-python-devcon AS mrd_converter

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Set the timezone to ensure good timestamps in logfiles
ENV TZ=America/Toronto

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git make curl sudo xvfb bzip2 gcc g++ python3-pip \
                    libglib2.0-0 libgl1 libxrender1 libxkbcommon-x11-0 libdbus-1-3 cmake && \
    apt-get clean

# Clone the python-ismrmrd-server repository
RUN mkdir -p /opt/code/python-ismrmrd-server && \
    git clone https://github.com/shimming-toolbox/python-ismrmrd-server.git /opt/code/python-ismrmrd-server

# Clone the spinal-cord-toolbox repository and install it
RUN mkdir -p /opt/code/spinal-cord-toolbox && \
    git clone --branch 7.0 --depth 1 https://github.com/spinalcordtoolbox/spinalcordtoolbox.git /opt/code/spinal-cord-toolbox && \
    cd /opt/code/spinal-cord-toolbox && \
    ./install_sct -y

# Clone the shimming-toolbox repository and install it
RUN mkdir -p /opt/code/shimming-toolbox && \
    git clone --branch 1.2 --depth 1 https://github.com/shimming-toolbox/shimming-toolbox.git /opt/code/shimming-toolbox && \
    cd /opt/code/shimming-toolbox && \
    make install PLUGIN=false

SHELL ["/bin/bash", "-c"]
# Activate automatically the conda environment when the container starts
RUN echo 'source /root/shimming-toolbox/python/etc/profile.d/conda.sh' >> ~/.bashrc

# Add to path
ENV PATH /root/shimming-toolbox/python/bin:$PATH

# Install specific versions
RUN pip install numpy==1.26.4 h5py==3.14.0

# Clone the mrd2nii repository and install it
RUN mkdir -p /opt/code/mrd2nii && \
    git clone https://github.com/shimming-toolbox/mrd2nii.git /opt/code/mrd2nii && \
    cd /opt/code/mrd2nii && \
    pip install .

# Install bet2
RUN conda install -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/ -c conda-forge fsl-bet2

# Install prelude
RUN mkdir /opt/code/prelude && \
    curl -o /opt/code/prelude/prelude -JL https://github.com/shimming-toolbox/binaries/raw/master/prelude && \
    sudo install /opt/code/prelude/prelude /usr/local/bin

FROM mrd_converter AS python-mrd-devcontainer

# Add ST modules
COPY ./python_modules/st_masking.py    /opt/code/python-ismrmrd-server
COPY ./python_modules/st_masking.json  /opt/code/python-ismrmrd-server
# Set the st_masking module as the default module
# CMD ["/bin/bash", "-c", "source /root/shimming-toolbox/python/etc/profile.d/conda.sh && conda activate shim-dev && python /opt/code/python-ismrmrd-server/main.py -v -H=0.0.0.0 -p=9002 -l=/tmp/python-ismrmrd-server.log -d=st_masking -s"]
CMD /bin/bash

# Clean up to reduce image size
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    rm -rf /root/.cache/pip && \
    rm -rf /var/cache/apt/archives && \
    apt-get autoremove -y && \
    apt-get autoclean
