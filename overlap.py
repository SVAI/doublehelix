#!/usr/bin/env python

import vcf
import csv
import sys
import collections
import itertools

chrom_ord = dict(zip(map(str, range(1, 23)), range(1, 23)))
chrom_ord['X'] = 24

def merged_read():
    reader_G91 = vcf.Reader(open('Data/G91716.vcf'))
    reader_G97 = vcf.Reader(open('Data/G97552.vcf'))

    record_G91 = next(reader_G91)
    record_G97 = next(reader_G97)

    while True:
        try:
            location_G91 = chrom_ord[record_G91.CHROM], record_G91.POS
            location_G97 = chrom_ord[record_G97.CHROM], record_G97.POS
            if location_G91 == location_G97:
                yield (record_G91, record_G97)
                record_G91 = next(reader_G91)
                record_G97 = next(reader_G97)
            elif location_G91 < location_G97:
                record_G91 = next(reader_G91)
            else:
                record_G97 = next(reader_G97)
        except StopIteration:
            break

def all_samples():
    def common_records():
        for r91, r97 in merged_read():
            if (r91.FILTER == [] and r91.FILTER == r97.FILTER
                and r91.REF == r97.REF and r91.ALT == r97.ALT):
                yield r91
                yield r97
    for r in common_records():
        for s in r.samples:
            yield s

out_file = sys.stdout
dwriter = csv.DictWriter(out_file, ['CHROM', 'POS', 'sample', 'called', 'gt'])
dwriter.writeheader()
for i, sample in enumerate(all_samples()):
    if i >= 1e10:
        break
    drow = {'CHROM': sample.site.CHROM,
            'POS': sample.site.POS,
            'sample': sample.sample,
            'called': sample.called,
            'gt': sample['GT']}
    dwriter.writerow(drow)

