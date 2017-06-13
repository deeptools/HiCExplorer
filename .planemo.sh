#!/bin/bash

tmp_dir=`mktemp -d`
conda create -c bioconda -c conda-forge --name python2.7 python=2.7.13 -y
source activate python2.7
# conda install -c bioconda -c conda-forge planemo -y
pip install planemo
planemo_bin='which planemo'
python_bin='which python'

# echo $planemo_bin
# source deactivate
# echo $CONDA_PREFIX
planemo database_create galaxy
# curl https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -o miniconda.sh
# chmod +x miniconda.sh
# ./miniconda.sh -b -p $tmp_dir/conda
# echo $CONDA_PREFIX
planemo conda_init --conda_prefix $tmp_dir/conda
# export PATH=$tmp_dir/conda/bin:$PATH
# conda install -y -c bioconda samtools python=2.7.13 numpy scipy matplotlib=2.0.0 nose flake8 pytables biopython pysam pybigwig intervaltree future six pandas
source deactivate
export PATH="$PATH_WITHOUT_CONDA"
hash -r

# source activate hicexplorer_galaxy
# echo $CONDA_PREFIX
pip install .

# echo $CONDA_PREFIX
# Galaxy wrapper testing
$planemo_bin $python_bin test --install_galaxy --galaxy_branch release_17.01 --skip_venv --conda_prefix $tmp_dir/prefix --conda_exec $tmp_dir/conda/bin/conda --conda_dependency_resolution --postgres galaxy/wrapper
# planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
# /home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
# source deactivate