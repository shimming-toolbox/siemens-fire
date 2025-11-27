@echo off
setlocal enabledelayedexpansion

rem This script takes a Docker image and creates a chroot image (.img)
rem Note that this script also requires docker_tar_to_chroot.sh to be located in the same folder
 
rem Syntax: docker_to_chroot.bat kspacekelvin/fire-python fire-python-chroot.img optional_buffer_size_in_mb
 
if "%1"=="" GOTO wrongargnum

set "SCRIPT_DIR=%~dp0"
call :get_parent_directory REPO_DIR "%SCRIPT_DIR:~0,-1%"
set DOCKER_NAME=st_image
set CHROOT_FILE=%1
set EXPORT_FILE=docker-export.tar

if "%3"=="" (
  set BUFFER_SIZE=1000
) else (
  set BUFFER_SIZE=%3
)

if exist %REPO_DIR%%EXPORT_FILE% (
  echo Warning -- %REPO_DIR%%EXPORT_FILE% exists and will be overwritten!
  del %REPO_DIR%%EXPORT_FILE%
)

rem Create a Docker image from the Dockerfile
echo ----------------------------------------------------------------------
echo Creating Docker image %DOCKER_NAME% from Dockerfile
echo ----------------------------------------------------------------------

docker build --no-cache --platform linux/amd64 -t %DOCKER_NAME% -f %SCRIPT_DIR%Dockerfile %REPO_DIR%
rem The --platform linux/amd64 option ensures the Docker image is built for the x86_64 (amd64) architecture.
rem This matches the architecture required by the MARS system hardware.

rem Create a Docker container and export to a .tar file
echo ------------------------------------------------------------
echo Exporting Docker image %DOCKER_NAME%
echo ------------------------------------------------------------
 
docker create --name tmpimage %DOCKER_NAME%
docker export -o %REPO_DIR%%EXPORT_FILE% tmpimage
docker rm tmpimage
docker rmi %DOCKER_NAME%
docker builder prune -a -f

rem Run a privileged Docker to create the chroot file
docker run -it --rm          ^
           --privileged=true ^
           -v "%REPO_DIR:~0,-1%":/share  ^
           ubuntu            ^
           /bin/bash -c "sed -i -e 's/\r//g' /share/docker/docker_tar_to_chroot.sh && /share/docker/docker_tar_to_chroot.sh /share/%EXPORT_FILE% /share/%CHROOT_FILE%.img !BUFFER_SIZE!"
del %EXPORT_FILE%

tar -acvf %CHROOT_FILE%.zip %CHROOT_FILE%.img

goto eof

:get_parent_directory
	set "%~1=%~dp2"
	exit /b 0

:wrongargnum
	echo Syntax: docker_to_chroot.bat chroot_file_name optional_buffer_size_in_mb

:eof
