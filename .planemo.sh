#!/bin/bash

tmp_dir=`mktemp -d`
pip install planemo
planemo database_create galaxy
planemo conda_init --conda_prefix $tmp_dir/conda

planemo test --install_galaxy --galaxy_branch release_17.01 --skip_venv --postgres galaxy/wrapper 

