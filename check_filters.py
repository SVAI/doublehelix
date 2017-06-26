#!/usr/bin/env python

import vcf
import collections

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


cnt = collections.Counter()
total = 0
for record_no, (r91, r97) in enumerate(merged_read()):
    if record_no >= 1e10:
        break
    assert r91.POS == r97.POS
    r91_pass = r91.FILTER == []
    r97_pass = r97.FILTER == []
    bucket = 'r91_pass = %d ; r97_pass = %d' % (r91_pass, r97_pass)
    cnt[bucket] += 1
    total += 1
print(total)
for k, v in cnt.items():
    print('%s ; %10d' % (k, v))
