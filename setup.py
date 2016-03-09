#-*- coding: utf-8 -*-

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
        print "This does not appear to be a Git repository."
        return
    try:
        p = subprocess.Popen(["git", "describe",
                              "--tags", "--always"],
                             stdout=subprocess.PIPE)
    except EnvironmentError:
        print "unable to run git, leaving hicexplorer/_version.py alone"
        return
    stdout = p.communicate()[0]
    if p.returncode != 0:
        print "unable to run git, leaving hicexplorer/_version.py alone"
        return
    ver = stdout.strip()
    f = open("hicexplorer/_version.py", "w")
    f.write(VERSION_PY % ver)
    f.close()
    print "set hicexplorer/_version.py to '%s'" % ver


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


class install(_install):
    def run(self):
        _install.run(self)
        if os.environ.get('DEEP_TOOLS_NO_CONFIG', False):
            return
        self.config_file = self.install_platlib + \
            "/hicexplorer/config/hicexplorer.cfg"

       # check installation of several components
        samtools_installed = self.checkProgramIsInstalled(
            'samtools', 'view',
            'http://samtools.sourceforge.net/',
            'correctGCbias')

        bedGraphToBigWig_installed = self.checkProgramIsInstalled(
            'bedGraphToBigWig', '-h',
            'http://hgdownload.cse.ucsc.edu/admin/exe/',
            'bamCoverage, bamCompare, correctGCbias')

        bigwigInfo_installed = self.checkProgramIsInstalled(
            'bigWigInfo', '-h',
            'http://hgdownload.cse.ucsc.edu/admin/exe/',
            'bigwigCompare')

        if not samtools_installed or not bedGraphToBigWig_installed \
                or not bigwigInfo_installed:

            msg = "\n##########\nSome tools were not fund.\n"\
                "If you already have a copy of this programs installed\n"\
                "please be sure that they are found in your PATH or\n"\
                "that they referred in the configuration file of hicexplorer\n"\
                "located at:\n\n {}\n\n".format(self.config_file)
            sys.stderr.write(msg)

    def checkProgramIsInstalled(self, program, args, where_to_download,
                                affected_tools):
        try:
            _out = subprocess.Popen([program, args], stderr=subprocess.PIPE,
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


########

setup(
    name='HiCExplorer',
    version=get_version(),
    author='Fidel Ramirez',
    author_email='ramirez@ie-freiburg.mpg.de',
    packages=['hicexplorer'],
    scripts=['bin/findRestSite', 'bin/hicBuildMatrix', 'bin/hicCorrectMatrix',
             'bin/hicCorrelate', 'bin/hicFindEnrichedContacts', 'bin/hicFindTADs',
             'bin/hicMergeMatrixBins', 'bin/hicPlotMatrix', 'bin/hicComputeSaturation',
             'bin/hicPlotTADs', 'bin/hicSumMatrices', 'bin/hicExport'],
    include_package_data=True,
    package_data={'': ['config/hicexplorer.cfg']},
    #url='http://pypi.python.org/pypi/hicexplorer/',
    license='LICENSE.txt',
    description='Set of programms to process, analyze and visualize Hi-C data'
    'BAM format.',
    long_description=open('README.rst').read(),
    install_requires=[
        "numpy >= 1.8.1",
        "scipy >= 0.14.0",
        "matplotlib >= 1.3.1",
        "pysam >= 0.8.3",
        "bx-python >= 0.7.1",
        "biopython >= 1.65"
    ],
    cmdclass={'sdist': sdist, 'install': install}
)
