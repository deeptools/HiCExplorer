#!/bin/bash

blah=`mktemp -d`
/home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo database_create galaxy
/home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo conda_init --conda_prefix $blah/conda
export PATH=$blah/conda/bin:$PATH
conda create -y -c bioconda --name hicexplorer_galaxy samtools python=2.7.13 numpy scipy matplotlib=1.5.3 nose flake8 pytables biopython pysam pybigwig intervaltree future six
source activate hicexplorer_galaxy

pip install .

# Galaxy wrapper testing
/home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo test --skip_venv --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
