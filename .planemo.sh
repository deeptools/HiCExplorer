#!/bin/bash

tmp_dir=`mktemp -d`
# conda create -c bioconda -c conda-forge --name python2.7 python=2.7.13 -y 
# source activate python2.7
pip install planemo

planemo database_create galaxy
planemo conda_init --conda_prefix $tmp_dir/conda
planemo conda_install galaxy/wrapper --conda_prefix $tmp_dir/conda
python setup.py install
planemo test --install_galaxy --galaxy_branch release_17.01 --skip_venv --conda_prefix $tmp_dir/conda --conda_dependency_resolution --postgres galaxy/wrapper
# source deactivate