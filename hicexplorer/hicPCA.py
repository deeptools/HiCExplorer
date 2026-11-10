"""hicPCA, implemented in C++.

main() runs the C++ executable hicPCA and parse_arguments() is built from its
--help-json (hicexplorer._cpp). correlateEigenvectorWithGeneTrack and
correlateEigenvectorWithHistonMarkTrack keep the name, the signature and the
behaviour they had in HiCExplorer 3.7.6; their bodies call
hicexplorer._hicx.hicPCA (cpp/python/src/api/hicPCA.cpp).
"""
import warnings
import logging
log = logging.getLogger(__name__)

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

from hicexplorer._hicx import hicPCA as _impl


def parse_arguments(args=None):
    """The argparse parser of hicPCA, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPCA')


def correlateEigenvectorWithGeneTrack(pMatrix, pEigenvector, pGeneTrack):
    '''
    This function correlates the eigenvectors per chromosome with the gene
    density. If the correlation is negative, the eigenvector values are
    multiplied with -1.
    '''
    return _impl.correlateEigenvectorWithGeneTrack(pMatrix, pEigenvector, pGeneTrack)


def correlateEigenvectorWithHistonMarkTrack(pEigenvector, bwTrack, chromosome,
                                            start, end, pHistonMarkTrack,
                                            pHistonMarkType):
    """
        This function flip the signs only if both compartments exist for the
        given `chromosome`, otherwise it doesn't change the signs and
        the investigation should be done manually.
    """
    return _impl.correlateEigenvectorWithHistonMarkTrack(pEigenvector, bwTrack, chromosome, start, end,
                                                         pHistonMarkTrack, pHistonMarkType)


def main(args=None):
    """Run the C++ implementation of hicPCA with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPCA')(args)
