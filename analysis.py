import pandas as pd
import numpy as np
import re
#pd.set_option('display.max_columns', None)

df = pd.read_csv ('/Users/kehuang/Desktop/hackathon/G97552.csv', delimiter = '\t')
ef = pd.read_csv ('/Users/kehuang/Desktop/hackathon/G91716.csv', delimiter = '\t')

# OF_112015SJIA_2 GERMLINE
# OF_010116NF2_a TUMOR
# NF2_XY_s BROTHER

#df s3 OF_010116NF2_a TUMOR
#ef s2 OF_112015SJIA_2 GERMLINE

def low_quality(x):
    if -4.7847 <= x < -2.7677:
        return "FAIL"
    elif -4.7847 <= x < -2.7677:
        return "FAIL"
    elif x < -38557.2084:
        return "FAIL"
    elif -38557.2084 <= x < -4.7847:
        return "FAIL"
    elif -2.6522 <= x < -1.21:
        return "FAIL"
    elif -4.7453 <= x < -2.6522:
        return "FAIL"
    elif x < -45335.9305:
        return "FAIL"
    elif -45335.9305 <= x < -4.7453:
        return "FAIL"
    else:
        return "PASS"

df['FILTER'] = df['info.VQSLOD'].apply(low_quality)
ef['FILTER'] = ef['info.VQSLOD'].apply(low_quality)

# OF_112015SJIA_2 GERMLINE
# OF_010116NF2_a TUMOR
# NF2_XY_s BROTHER

#df s3 OF_010116NF2_a TUMOR
#ef s2 OF_112015SJIA_2 GERMLINE

tumor_col = [i for i in df.columns if not i.startswith('s1') and not i.startswith('s2')]
germline_col = [i for i in ef.columns if not i.startswith('s1') and not i.startswith('s3')]

tumor = df[df.FILTER == 'PASS'][tumor_col]
germline = ef[ef.FILTER == 'PASS'][germline_col]

of = pd.merge(tumor, germline, how='outer', on=['CHROM', 'POS'], suffixes = ('_t', '_g'))

of[(of['REF_g'] == of['REF_t']) & (of['ALT_g'] == of['ALT_t']) \
   & (~of['s2'].isnull()) & (~of['s3'].isnull()) \
   & (of['s3.GT'] <> of['s2.GT']) \
  ][['CHROM', 'POS', 's2', 's3', 'REF_g', 'ALT_g', 's2.GT', 's3.GT']].to_csv('/Users/kehuang/Desktop/hackathon/diff.csv')
