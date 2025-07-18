FROM kspacekelvin/fire-python-devcon

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git make curl sudo xvfb bzip2 gcc g++ python3-pip \
                    libglib2.0-0 libgl1 libxrender1 libxkbcommon-x11-0 libdbus-1-3 && \
    apt-get clean

# Clone the python-ismrmrd-server repository
RUN mkdir -p /opt/code/python-ismrmrd-server && \
    git clone https://github.com/kspaceKelvin/python-ismrmrd-server.git /opt/code/python-ismrmrd-server
# Remove unecessary files
RUN rm -rf /opt/code/python-ismrmrd-server/.git && \
    rm -rf /opt/code/python-ismrmrd-server/.github/workflows && \
    rm -rf /opt/code/python-ismrmrd-server/.gitattributes && \
    rm -rf /opt/code/python-ismrmrd-server/.gitignore && \
    rm -rf /opt/code/python-ismrmrd-server/.vscode && \
    rm -rf /opt/code/python-ismrmrd-server/custom && \
    rm -rf /opt/code/python-ismrmrd-server/doc && \
    rm -rf /opt/code/python-ismrmrd-server/docker && \
    rm -rf /opt/code/python-ismrmrd-server/RunClientServerRecon.ipynb && \
    rm -rf /opt/code/python-ismrmrd-server/RunServerRecon.ipynb && \
    rm -rf /opt/code/python-ismrmrd-server/mrd2gif.py && \
    rm -rf /opt/code/python-ismrmrd-server/invertcontrast.py && \
    rm -rf /opt/code/python-ismrmrd-server/invertcontrast.json && \
    rm -rf /opt/code/python-ismrmrd-server/simplefft.py && \
    rm -rf /opt/code/python-ismrmrd-server/analyzeflow.py && \
    rm -rf /opt/code/python-ismrmrd-server/report.py

# Clone the spinal-cord-toolbox repository and install it
RUN mkdir -p /opt/code/spinal-cord-toolbox && \
    git clone --branch 7.0 --depth 1 https://github.com/spinalcordtoolbox/spinalcordtoolbox.git /opt/code/spinal-cord-toolbox && \
    cd /opt/code/spinal-cord-toolbox && \
    ./install_sct -y
# Remove unecessary files
RUN rm -rf /opt/code/spinal-cord-toolbox/.git && \
    rm -rf /opt/code/spinal-cord-toolbox/.github && \
    rm -rf /opt/code/spinal-cord-toolbox/.gitignore && \
    rm -rf /opt/code/spinal-cord-toolbox/.ci.sh && \
    rm -rf /opt/code/spinal-cord-toolbox/.readthedocs.yaml && \
    rm -rf /opt/code/spinal-cord-toolbox/contrib && \
    rm -rf /opt/code/spinal-cord-toolbox/documentation && \
    rm -rf /opt/code/spinal-cord-toolbox/testing && \
    rm -rf /opt/code/spinal-cord-toolbox/CONTRIBUTING.rst && \
    rm -rf /opt/code/spinal-cord-toolbox/MANIFEST.in && \
    rm -rf /opt/code/spinal-cord-toolbox/batch_processing.sh && \
    rm -rf /opt/code/spinal-cord-toolbox/install_sct && \
    rm -rf /opt/code/spinal-cord-toolbox/install_sct.bat && \
    rm -rf /opt/code/spinal-cord-toolbox/requirements.txt && \
    rm -rf /opt/code/spinal-cord-toolbox/setup.cfg && \
    rm -rf /opt/code/spinal-cord-toolbox/setup.py

# Clone the shimming-toolbox repository and install it
RUN mkdir -p /opt/code/shimming-toolbox && \
    git clone --branch 1.2 --depth 1 https://github.com/shimming-toolbox/shimming-toolbox.git /opt/code/shimming-toolbox && \
    cd /opt/code/shimming-toolbox && \
    make install PLUGIN=false
# Create the conda environment
SHELL ["/bin/bash", "-c"]
RUN cd /opt/code/shimming-toolbox/shimming-toolbox && \
    source /root/shimming-toolbox/python/etc/profile.d/conda.sh && \
    conda create -n shim-dev -c conda-forge dcm2niix python=3.10 && \
    conda activate shim-dev && \
    pip install numpy==1.26.4 h5py==3.14.0 && \
    pip install -e ."[docs,dev]"
# Activate automatically the conda environment when the container starts
RUN echo 'source root/shimming-toolbox/python/etc/profile.d/conda.sh' >> ~/.bashrc && \
    echo 'conda activate shim-dev' >> ~/.bashrc
# Remove unecessary files
RUN rm -rf /opt/code/shimming-toolbox/.git && \
    rm -rf /opt/code/shimming-toolbox/.github && \
    rm -rf /opt/code/shimming-toolbox/.gitignore && \
    rm -rf /opt/code/shimming-toolbox/.coveragerc && \
    rm -rf /opt/code/shimming-toolbox/.pre-commit-config.yaml && \
    rm -rf /opt/code/shimming-toolbox/docs && \
    rm -rf /opt/code/shimming-toolbox/examples && \
    rm -rf /opt/code/shimming-toolbox/installer && \
    rm -rf /opt/code/shimming-toolbox/Makefile && \
    rm -rf /opt/code/shimming-toolbox/readthedocs.yml && \
    rm rm -rf /opt/code/shimming-toolbox/shimming-toolbox/test && \
    rm -rf /opt/code/shimming-toolbox/shimming-toolbox/.github && \
    rm -rf /opt/code/shimming-toolbox/shimming-toolbox/MANIFEST.in && \
    rm -rf /opt/code/shimming-toolbox/shimming-toolbox/conftest.py && \
    rm -rf /opt/code/shimming-toolbox/shimming-toolbox/pytest.ini && \
    rm -rf /opt/code/shimming-toolbox/shimming-toolbox/requirements_st.txt && \
    rm -rf /opt/code/shimming-toolbox/shimming-toolbox/setup.py

# Clone the mrd2nii repository and install it
RUN mkdir -p /opt/code/mrd2nii && \
    git clone https://github.com/shimming-toolbox/mrd2nii.git /opt/code/mrd2nii && \
    cd /opt/code/mrd2nii && \
    source /root/shimming-toolbox/python/etc/profile.d/conda.sh && \
    source /root/shimming-toolbox/python/bin/activate && \
    conda activate shim-dev && \
    pip install -e .

# Install bet2
RUN source /root/shimming-toolbox/python/etc/profile.d/conda.sh && \
    source /root/shimming-toolbox/python/bin/activate && \
    conda activate shim-dev && \
    conda install -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/ -c conda-forge fsl-bet2
    
# Install prelude
RUN mkdir /opt/code/prelude && \
    curl -o /opt/code/prelude/prelude -JL https://github.com/shimming-toolbox/binaries/raw/master/prelude && \
    sudo install /opt/code/prelude/prelude /usr/local/bin

# Add ST modules
COPY ./python_modules/st_masking.py    /opt/code/python-ismrmrd-server
COPY ./python_modules/st_masking.json  /opt/code/python-ismrmrd-server
# Set the st_masking module as the default module
CMD [ "python3", "/opt/code/python-ismrmrd-server/main.py", "-v", "-H=0.0.0.0", "-p=9002", "-l=/tmp/python-ismrmrd-server.log", "-d=st_masking"]

# Clean up to reduce image size
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    rm -rf /root/.cache/pip && \
    rm -rf /var/cache/apt/archives && \
    apt-get autoremove -y && \
    apt-get autoclean
