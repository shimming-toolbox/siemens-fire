# Siemens-FIRE

This repository contains tools, scripts, and instructions for building a chroot disk image of _Shimming Toolbox_ and deploying it on Siemens MARS MRI systems using the _Framework for Image REconstruction_ (FIRE). This approach simplifies advanced shimming protocols by removing the need for an external laptop running _Shimming Toolbox_. It provides a step-by-step guide to create the image locally with Docker, compress it, and transfer it to a USB drive for installation on the scanner hardware. The repository also includes guidelines for testing Python modules both locally and directly on the MARS system.

Many of the scripts are adapted from [Kelvin Chow's _python-ismrmrd-server_ repository](https://github.com/kspaceKelvin/python-ismrmrd-server). Please refer to that repository’s documentation for additional details.

---

## Table of Contents
* 1. [Creating a chroot disk image for Shimming Toolbox](#creating-a-chroot-disk-image-for-shimming-toolbox)
  * 1.1 [Prerequisites](#prerequisites)
  * 1.2 [Overview](#overview)
  * 1.3 [Step 1 - Clone and navigate to the repository](#step-1---clone-and-navigate-to-the-repository)
  * 1.4 [Step 2 - Create the chroot image](#step-2---create-the-chroot-image)
  * 1.5 [Step 3 - Export the chroot image to your USB key](#step-3---export-the-chroot-image-to-your-usb-key)
  * 1.6 [File Size Reference](#file-size-reference)
  * 1.7 [Next Steps](#next-steps)
* 2. [Python modules](#python-modules)
* 3. [Testing a python module locally](#3-testing-a-python-module-locally)
* 4. [Testing a python module on MARS system hardware](#4-testing-a-python-module-on-mars-system-hardware)

---

<a name="creating-a-chroot-disk-image-for-shimming-toolbox"></a>
## 1. Creating a chroot disk image for Shimming Toolbox
This guide walks you through creating on your computer a disk image (.img) containing a complete Linux filesystem intended to be used as the root directory for a changed root (chroot) environment. The resulting image can be exported and deployed on MARS system hardware to run Shimming Toolbox locally on the MRI scanner.

<a name="prerequisites"></a>
### 1.1 Prerequisites
Before starting, ensure you have:
- **Docker**: Installed and running on your system
- **Disk space**: Between 30-35GB free disk space for the build process
- **USB drive**: At least 8GB capacity for image transfer
- **Time**: Allow 30-40 minutes for the complete process
- **Permissions**: Root/sudo access may be required for some operations

<a name="overview"></a>
### 1.2 Overview
The process consists of three main steps:
1. Build a Docker image containing Shimming Toolbox
2. Export the container filesystem and create a chroot-compatible disk image
3. Compress and transfer the image to a USB drive for deployment

<a name="step-1---clone-and-navigate-to-the-repository"></a>
### 1.3 Step 1 - Clone and navigate to the repository
```bash
git clone https://github.com/shimming-toolbox/siemens-fire.git <path/to/repository>
cd <path/to/repository>
```
Verify you're in the correct directory by checking for the required files:
```bash
ls -la
```
You should see `Dockerfile` and the build scripts (`docker_to_chroot.sh`, `tar_to_chroot.sh`).

<a name="step-2---create-the-chroot-image"></a>
### 1.4 Step 2 - Create the chroot image
Make the build scripts executable and run the image creation process:
```bash
chmod +x ./docker_to_chroot.sh
chmod +x ./tar_to_chroot.sh
./docker_to_chroot.sh st_image st_chroot.img
```
**What happens during this step:**
1. **Docker build**: Creates a Docker image containing Shimming Toolbox and its dependencies
2. **Container export**: Exports the container's filesystem to a tarball
3. **Image creation**: Converts the tarball to a Linux root filesystem (.img)
> [!NOTE]
> **Expected duration:** 14 minutes
> **Output file:** `st_chroot.img` (approximately 8.8GB)

<a name="step-3---export-the-chroot-image-to-your-usb-key"></a>
### 1.5 Step 3 - Export the chroot image to your USB key
**Step 3.1 - Prepare and verify**

Connect your USB drive and identify its mount point:
```bash
# On macOS
df -h | grep /Volumes/

# On Linux
df -h | grep /media/
```
**Step 3.2 - Create compressed archive**

Compress the image to reduce transfer time and storage requirements:
```bash
zip st_chroot.zip st_chroot.img
```
> [!NOTE]
> **Output file:** `st_chroot.zip` (approximately 2.6GB)

**Step 3.3 - Transfer to your USB drive**

Replace `/Volumes/YOUR_USB_NAME/` with your actual USB drive path:
```bash
cp st_chroot.zip </Volumes/YOUR_USB_NAME/>
```
> [!NOTE]
> **Expected transfer time:** 2-10 minutes depending on USB drive speed (USB 3.0+ recommended)

**Step 3.4 - Verify and cleanup**

Verify the transfer completed successfully:
```bash
ls -la </Volumes/YOUR_USB_NAME>/st_chroot.zip
```
Compare file sizes to ensure integrity:
```bash
ls -l st_chroot.zip </Volumes/YOUR_USB_NAME>/st_chroot.zip
```
Only after successful verification, clean up local files:
```bash
rm st_chroot.zip
rm st_chroot.img
```
> [!WARNING]
> Don't delete the local files until you've verified the USB transfer completed successfully.

<a name="file-size-reference"></a>
### 1.6 File Size Reference
| File | Approximate Size |
|------|------------------|
| st_chroot.img | 8-9 GB |
| st_chroot.zip | 2-3 GB |
| Required USB space | 8+ GB |

<a name="next-steps"></a>
### 1.7 Next Steps
After completing this process:
1. Safely eject the USB drive
2. Transport to the MRI scanner system
3. Follow the deployment instructions for mounting the chroot environment on the MARS hardware

---

<a name="python-modules"></a>
## 2. Python modules
To perform shimming with _Shimmin Toolbox_ on the scanner, different python modules are linked to each acquisition : 

1. Localizer : No python module
2. MP-RAGE : Creation of the masked ROI with `st_masking.py`
3. GRE : Creation of the fieldmap with `st_fieldmapping.py`
4. Dummy EPI : Compute the currents with `st_optimizing.py`
5. Shimmed EPI : Shim with `st_shimming.py`

Each of these python modules take the raw data from the acquisition, converts it to a NIfTI file, runs 
the appropriate _Shimming Toolbox_ command, then converts it back to raw data that is sent back to the scanner.
Before an acquisition, to change de parameters of a python module, please update the corresponding JSON sidecar (e.g.
`st_shimming.json` for `st_shimming.py`).

---

<a name="3-testing-a-python-module-locally"></a>
## 3. Testing a python module locally
This guide describes how to build, launch, and interact with a Docker image for processing .h5 datasets on your computer using a Siemens-fire server and client CLI.

1. Navigate to the _siemens-fire_ directory
```bash
cd <path/to/repository>
```

2. Build the Docker image
```bash
docker build --platform linux/amd64 -t test_st .
```

3. Start the server on a new terminal window
```bash
docker run -p=9002:9002 --rm -it -v /tmp:/tmp test_st 
```
> [!NOTE]
> The server can be stopped any time using `control + c`

5. Start an interactive Bash CLI on a new terminal window
```bash
docker run --rm -it -v /tmp:/tmp test_st /bin/bash
```
> [!NOTE]
> The interactive CLI can be stopped by typing `exit`

6. Ensure you have test .h5 data files in your computer's /tmp directory
7. In the interactive CLI, launch the client
```bash
python3 /opt/code/python-ismrmrd-server/client.py -a host.docker.internal -c <module> -p 9002 -G "dataset" -o <path/to/output> <path/to/input>
```

For more details, please consult [Kelvin Chow's documentation](https://github.com/kspaceKelvin/python-ismrmrd-server/blob/master/doc/docker.md) on this subject.

---

<a name="4-testing-a-python-module-on-mars-system-hardware"></a>
## 4. Testing a python module on MARS system hardware
> [!NOTE]
> **TODO**

