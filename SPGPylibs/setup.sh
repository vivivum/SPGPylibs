#!/bin/bash

#compile the rte code

cd ./PHItools/cmilos/lib
make clear
make

cd ..
make clear
make

cd ../../

#conda env create -f environment.yml
#pip install -r requrements.txt

source activate hrt_pipeline_env

