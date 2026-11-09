#!/usr/bin/env python3
"""Writes the BAM pair of the hicBuildMatrix SA tag test (issue #387).

    python cpp/scripts/make_sa_semicolon_bams.py hicexplorer/test/test_data/hicBuildMatrix/sa_semicolon

The reference 'Sc1;HRSCAF=2' contains a ';'. The first read of the first pair of
R1 carries an SA tag "Sc1;HRSCAF=2,3000,+,50S50M,60,0;" that names one
supplementary alignment, and that record follows it in the file. Counting the
';' characters of the tag gives 2 and makes the reader consume the next pair.
Also writes the restriction site file.
"""
import os
import sys

import pysam

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
names = ['Sc1;HRSCAF=2', 'chr2']
header = {'HD': {'VN': '1.4', 'SO': 'unsorted'},
          'SQ': [{'SN': n, 'LN': 100000} for n in names]}


def record(bam, qname, flag, ref, pos, mapq=60, cigar='100M', sa=None):
    read = pysam.AlignedSegment(bam.header)
    read.query_name = qname
    read.flag = flag
    read.reference_id = names.index(ref)
    read.reference_start = pos
    read.mapping_quality = mapq
    read.cigarstring = cigar
    qlen = sum(int(c) for c in cigar.replace('S', 'M').split('M')[:-1])
    read.query_sequence = ('ACGTTGCA' * 20)[:qlen]
    read.query_qualities = pysam.qualitystring_to_array('I' * qlen)
    if sa:
        read.set_tag('SA', sa)
    return read


pairs = [
    ('p1', 'Sc1;HRSCAF=2', 1000, 'chr2', 1000),
    ('p2', 'chr2', 5000, 'Sc1;HRSCAF=2', 8000),
    ('p3', 'Sc1;HRSCAF=2', 20000, 'chr2', 30000),
    ('p4', 'chr2', 40000, 'chr2', 60000),
    ('p5', 'Sc1;HRSCAF=2', 50000, 'Sc1;HRSCAF=2', 70000),
]
with pysam.AlignmentFile(os.path.join(out, 'R1.bam'), 'wb', header=header) as r1, \
        pysam.AlignmentFile(os.path.join(out, 'R2.bam'), 'wb', header=header) as r2:
    for name, ref1, pos1, ref2, pos2 in pairs:
        sa = 'Sc1;HRSCAF=2,3000,+,50S50M,60,0;' if name == 'p1' else None
        r1.write(record(r1, name, 65, ref1, pos1, sa=sa))
        if sa:
            r1.write(record(r1, name, 65 + 2048, 'Sc1;HRSCAF=2', 3000, cigar='50S50M'))
        r2.write(record(r2, name, 129, ref2, pos2))

with open(os.path.join(out, 'sites.bed'), 'w') as handle:
    for n in names:
        for start in range(0, 100000, 500):
            handle.write('{}\t{}\t{}\t.\t0\t+\n'.format(n, start + 250, start + 256))
