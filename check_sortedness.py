#!/usr/bin/env python

import sys
import vcf

chrom_ord = dict(zip(map(str, range(1, 23)), range(1, 23)))
chrom_ord['X'] = 24

reader = vcf.Reader(open(sys.argv[1]))
prev_record = next(reader)
record_no = 1
for record_no, record in enumerate(reader, 2):
    location = chrom_ord[record.CHROM], record.POS
    prev_location = chrom_ord[prev_record.CHROM], prev_record.POS
    if location <= prev_location:
        raise ValueError
    prev_record = record
print('All OK')
