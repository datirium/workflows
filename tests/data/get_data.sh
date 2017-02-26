#!/bin/bash
echo 'Get data from https://github.com/michael-kotliar/workflows_test_data/archive/master.zip'
wget https://github.com/michael-kotliar/workflows_test_data/archive/master.zip
echo 'Extract data'
unzip master.zip
echo 'Delete unused files'
mv workflows_test_data-master/* .
rm -rf workflows_test_data-master
rm master.zip
