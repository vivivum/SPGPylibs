#!/bin/bash

#RTE compilation
cd ./SPGPylibs/PHItools/cmilos/lib
make clear
#by default it compiles cmilos and py milos. Remove cython-build in case you do not want the py milos version
make 
make cython-build

cd ../../../

#conda env create -f environment.yml 
#pip install -r requirements.txt
#source activate SPGPylibs_env
