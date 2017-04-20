# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import re

from setuptools import setup
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
        update_version_py()
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)

# Install class to check for external dependencies from OS environment


class install(_install):

    def run(self):
        update_version_py()
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


setup(
    name='HiCExplorer',
    version=get_version(),
    author='Fidel Ramirez, Vivek Bhardwaj, Björn Grüning',
    author_email='deeptools@googlegroups.com',
    packages=['hicexplorer'],
    scripts=['bin/findRestSite', 'bin/hicBuildMatrix', 'bin/hicCorrectMatrix',
             'bin/hicCorrelate', 'bin/hicFindEnrichedContacts', 'bin/hicFindTADs',
             'bin/hicMergeMatrixBins', 'bin/hicPlotMatrix', 'bin/hicPlotDistVsCounts',
             'bin/hicPlotTADs', 'bin/hicSumMatrices', 'bin/hicExport', 'bin/hicInfo', 'bin/hicexplorer',
             'bin/hicQC', 'bin/hicCompareMatrices'],
    include_package_data=True,
    package_dir={'hicexplorer': 'hicexplorer'},
    package_data={'hicexplorer': ['qc_template.html']},
    url='http://hicexplorer.readthedocs.io',
    license='LICENSE.txt',
    description='Set of programs to process, analyze and visualize Hi-C data',
    long_description=open('README.rst').read(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=[
        "numpy >= 1.10.4",
        "scipy >= 0.17.1",
        "matplotlib >= 1.5.3",
        "pysam >= 0.8.3",
        "intervaltree >= 2.1.0",
        "biopython >= 1.65",
        "tables >= 3.2.2",
        "pandas >= 0.19.2",
        "pyBigWig >=0.2.8"],
    zip_safe=False,
    cmdclass={'sdist': sdist, 'install': install}
)
