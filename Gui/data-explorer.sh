# Startup script for running the packaged Data Explorer on Linux
cd bin
LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH ./ogs-gui
