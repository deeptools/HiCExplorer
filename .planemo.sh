#!/bin/bash

blah=`mktemp -d`
/home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo database_create galaxy
/home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo conda_init --conda_prefix $blah/conda
export PATH=$blah/conda/bin:$PATH
conda create -y -c bioconda --name deeptools_galaxy numpy scipy matplotlib=1.5.3 nose flake8 pytables biopython pysam pyBigWig intervaltree
source activate deeptools_galaxy
conda install -c bioconda samtools
#git clone --depth 1 https://github.com/galaxyproject/galaxy.git clone
#cd clone
#Add the custom data types
#sed -i '4i\    <datatype extension="deeptools_compute_matrix_archive" type="galaxy.datatypes.binary:CompressedArchive" subclass="True" display_in_upload="True"/>' config/datatypes_conf.xml.sample
#sed -i '5i\    <datatype extension="deeptools_coverage_matrix" type="galaxy.datatypes.binary:CompressedArchive" subclass="True" display_in_upload="True"/>' config/datatypes_conf.xml.sample
#./scripts/common_startup.sh --skip-venv --dev-wheels
#cd ..
#conda uninstall -y sqlite
pip install . 
/home/travis/build/maxplanck-ie/HiCExplorer/foo/bin/planemo test --install_galaxy --no_conda_auto_install --no_conda_auto_init --galaxy_branch release_17.01 --postgres galaxy/wrapper/
