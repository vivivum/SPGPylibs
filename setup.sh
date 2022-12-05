#!/bin/bash

#compile the rte code

cd ./src/SPGPylibs/PHItools/cmilos
make clean
make

cd -

#conda env create -f environment.yml
#pip install -r requrements.txt
#source activate hrt_pipeline_env

