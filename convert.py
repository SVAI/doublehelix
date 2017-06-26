import csv
import vcf

file = 'G97552'

with open('/Users/kehuang/Desktop/hackathon/' + file + '.csv', 'wb') as csvfile:

    out = csv.writer(csvfile, delimiter='\t')
    reader = vcf.Reader(open('/Users/kehuang/Desktop/hackathon/' + file + '.vcf', 'r'))

    formats = reader.formats.keys()
    infos = reader.infos.keys()
    other = ['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL'] 


    header = ['info.' + x for x in infos] + other \
            + ['s1'] + ['s1.' + x for x in formats] \
            + ['s2'] + ['s2.' + x for x in formats] \
            + ['s3'] + ['s3.' + x for x in formats]


    out.writerow(header)


    def flatten(x):
        if type(x) == type([]):
            x = ','.join(map(str, x))
        return x

    for record in reader:
        row = [flatten(record.INFO.get(x, None)) for x in infos]
        row += [record.CHROM, record.POS, record.REF, record.ALT, record.ID, record.QUAL]

        for sample in record.samples:
            row += [sample.sample]
            row += [flatten(getattr(sample.data, x, None)) for x in formats]

        out.writerow(row)

