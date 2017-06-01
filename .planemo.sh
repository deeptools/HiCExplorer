#!/bin/bash

tmp_dir=`mktemp -d`

planemo_test_env/bin/planemo database_create galaxy
planemo_test_env/bin/planemo conda_init --conda_prefix $tmp_dir/conda
export PATH=$tmp_dir/conda/bin:$PATH
conda create -y -c bioconda --name hicexplorer_galaxy samtools python=2.7.13 numpy scipy matplotlib=1.5.3 nose flake8 pytables biopython pysam pybigwig intervaltree
source activate hicexplorer_galaxy

pip install .

ls -lah
ls -lah planemo_test_env
ls -lah planemo_test_env/bin

# Galaxy wrapper testing
planemo_test_env/bin/planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
# /home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
