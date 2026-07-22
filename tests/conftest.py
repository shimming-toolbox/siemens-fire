import shutil
import subprocess
import time
import socket
import pytest
import os
import time

from python_modules import __PATH_TESTING_DATA__, __TMP_SHARE_SAVEDDATA__, __TMP_SHARE_DEBUG__

SERVER_PORT = 9020
SERVER_URL = f"http://localhost:{SERVER_PORT}"


"""
sudo lsof -i :9020
kill -9 PID
"""

@pytest.fixture(scope="function", autouse=True)
def clean_savedata_folder():
    """Cleans the savedata folder before each test to ensure a clean state."""
    if os.path.exists(__TMP_SHARE_SAVEDDATA__):
        shutil.rmtree(__TMP_SHARE_SAVEDDATA__)
    os.makedirs(__TMP_SHARE_SAVEDDATA__, exist_ok=False)


@pytest.fixture(scope="function", autouse=True)
def clean_debug_folder():
    """Cleans the debug folder before each test to ensure a clean state."""
    if os.path.exists(__TMP_SHARE_DEBUG__):
        shutil.rmtree(__TMP_SHARE_DEBUG__)
    os.makedirs(__TMP_SHARE_DEBUG__, exist_ok=False)


@pytest.fixture(scope="session", autouse=True)
def run_server():
    """Starts the application server process in the background for testing."""

    if not os.path.exists(os.path.join(__PATH_TESTING_DATA__, "log")):
        os.makedirs(os.path.join(__PATH_TESTING_DATA__, "log"))

    fname_log_file = f"{__PATH_TESTING_DATA__}/log/server.log"
    if os.path.exists(fname_log_file):
        os.remove(fname_log_file)

    process = subprocess.Popen(
        ["/root/shimming-toolbox/python/bin/python",
         "/opt/code/python-ismrmrd-server/main.py",
         "-H", "0.0.0.0",
         "-p", f"{SERVER_PORT}",
         "-l", fname_log_file,
         # "-v",
         "-s"],
        text=True
    )

    # Better way to do this?
    time.sleep(1)  # Give the server a moment to start up
        
    yield  # Hand control back to client test files
    
    # Session teardown: Terminate the background process gracefully
    process.terminate()
    process.wait()
