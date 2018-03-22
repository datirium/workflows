#!/usr/bin/bash

sed -i 's/cwd = os.getcwd()/cwd = os.path.dirname(__file__)/g' /opt/rose/ROSE_main.py
sed -i 's/cwd = os.getcwd()/cwd = os.path.dirname(__file__)/g' /opt/rose/ROSE_main_turbo.py
sed -i 's/ROSE_bamToGFF.py/%s\/ROSE_bamToGFF.py/g' /opt/rose/ROSE_main.py
sed -i 's/(nBin/(os.path.dirname(__file__),nBin/g' /opt/rose/ROSE_main.py
sed -i 's/ROSE_bamToGFF.py/%s\/ROSE_bamToGFF.py/g' /opt/rose/ROSE_main_turbo.py
sed -i 's/(nBin/(os.path.dirname(__file__),nBin/g' /opt/rose/ROSE_main_turbo.py
sed -i 's/ROSE_callSuper.R/%s\/ROSE_callSuper.R/g' /opt/rose/ROSE_main.py
sed -i 's/ROSE_callSuper.R/%s\/ROSE_callSuper.R/g' /opt/rose/ROSE_main_turbo.py
sed -i 's/controlName)/controlName,os.path.dirname(__file__))/g' /opt/rose/ROSE_main.py
sed -i 's/controlName)/controlName,os.path.dirname(__file__))/g' /opt/rose/ROSE_main_turbo.py

