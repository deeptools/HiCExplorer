# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import re

from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.
__version__ = '%s'
"""


def update_version_py():
    if not os.path.isdir(".git"):
        print("This does not appear to be a Git repository.")
        return
    try:
        p = subprocess.Popen(["git", "describe",
                              "--tags", "--always"],
                             stdout=subprocess.PIPE)
    except EnvironmentError:
        print("unable to run git, leaving hicexplorer/_version.py alone")
        return
    stdout = p.communicate()[0]
    if p.returncode != 0:
        print("unable to run git, leaving hicexplorer/_version.py alone")
        return
    ver = stdout.strip()
    f = open("hicexplorer/_version.py", "w")
    f.write(VERSION_PY % ver)
    f.close()
    print("set hicexplorer/_version.py to '%s'" % ver)


def get_version():
    try:
        f = open("hicexplorer/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


class sdist(_sdist):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)

# Install class to check for external dependencies from OS environment


class install(_install):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return

    def checkProgramIsInstalled(self, program, args, where_to_download,
                                affected_tools):
        try:
            subprocess.Popen([program, args],
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE)
            return True
        except EnvironmentError:
            # handle file not found error.
            # the config file is installed in:
            msg = "\n**{0} not found. This " \
                  "program is needed for the following "\
                  "tools to work properly:\n"\
                  " {1}\n"\
                  "{0} can be downloaded from here:\n " \
                  " {2}\n".format(program, affected_tools,
                                  where_to_download)
            sys.stderr.write(msg)

        except Exception as e:
            sys.stderr.write("Error: {}".format(e))


install_requires_py = ["numpy >= 1.19",
                       "scipy >= 1.10",
                       "matplotlib >= 3.6",
                       "ipykernel >= 6.25.2",
                       "pysam >= 0.21",
                       "intervaltree >= 3.1",
                       "biopython",
                       "tables >= 3.8",
                       "pandas >= 2.0",
                       "pybigwig >= 0.3",
                       "jinja2 >= 3.1.2",
                       "unidecode >= 1.3",
                       "hicmatrix >= 17",
                       "hic2cool >= 0.8.3",
                       "psutil >= 5.9",
                       "pygenometracks >= 3.8",
                       "fit_nbinom >= 1.2",
                       "cooler >= 0.9.3",
                       "krbalancing >= 0.0.5",
                       "pybedtools >= 0.9",
                       "future >= 0.18",
                       "tqdm >= 4.66",
                       "hyperopt >= 0.2.7",
                       "graphviz >= 0.20",
                       "scikit-learn >= 1.3,<1.4",
                       "imbalanced-learn >= 0.11",
                       "cleanlab >= 2.5"
                       ]


setup(
    name='HiCExplorer',
    version=get_version(),
    author='Joachim Wolff, Leily Rabbani, Bjoern Gruening, Vivek Bhardwaj, Fidel Ramírez',
    author_email='deeptools@googlegroups.com',
    packages=find_packages(),
    scripts=['bin/hicFindRestSite', 'bin/hicAggregateContacts', 'bin/hicBuildMatrix', 'bin/hicCorrectMatrix',
             'bin/hicCorrelate', 'bin/hicFindTADs', 'bin/hicMergeMatrixBins', 'bin/hicPlotMatrix', 'bin/hicPlotDistVsCounts',
             'bin/hicPlotTADs', 'bin/hicSumMatrices', 'bin/hicInfo', 'bin/hicexplorer',
             'bin/hicQC', 'bin/hicCompareMatrices', 'bin/hicPCA', 'bin/hicTransform', 'bin/hicPlotViewpoint',
             'bin/chicViewpointBackgroundModel', 'bin/chicPlotViewpoint', 'bin/chicViewpoint',
             'bin/chicAggregateStatistic', 'bin/chicDifferentialTest', 'bin/chicQualityControl', 'bin/chicSignificantInteractions',
             'bin/hicConvertFormat', 'bin/hicAdjustMatrix', 'bin/hicNormalize',
             'bin/hicAverageRegions', 'bin/hicPlotAverageRegions', 'bin/hicDetectLoops', 'bin/hicValidateLocations', 'bin/hicMergeLoops',
             'bin/hicCompartmentalization', 'bin/hicQuickQC', 'bin/hicPlotSVL', 'bin/hicCreateThresholdFile', 'bin/hicHyperoptDetectLoops',
             'bin/hicHyperoptDetectLoopsHiCCUPS', 'bin/hicMergeDomains', 'bin/hicDifferentialTAD', 'bin/chicExportData', 'bin/hicInterIntraTAD',
             'bin/hicTADClassifier', 'bin/hicTrainTADClassifier'
             ],
    include_package_data=True,
    package_dir={'hicexplorer': 'hicexplorer'},
    package_data={'hicexplorer': ['qc_template.html']},
    url='http://hicexplorer.readthedocs.io',
    license='LICENSE',
    description='Set of programs to process, analyze and visualize Hi-C data',
    long_description=open('README.rst').read(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=install_requires_py,
    zip_safe=False,
    cmdclass={'sdist': sdist, 'install': install}
)
