#!/bin/bash

tmp_dir=`mktemp -d`
source activate python2.7
planemo_bin='which planemo'
source deactivate
$planemo_bin database_create galaxy
$planemo_bin conda_init --conda_prefix $tmp_dir/conda
export PATH=$tmp_dir/conda/bin:$PATH
conda create -y -c bioconda --name hicexplorer_galaxy samtools python=2.7.13 numpy scipy matplotlib=1.5.3 nose flake8 pytables biopython pysam pybigwig intervaltree future six pandas

source activate hicexplorer_galaxy

pip install .


# Galaxy wrapper testing
$planemo_bin test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
# /home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
