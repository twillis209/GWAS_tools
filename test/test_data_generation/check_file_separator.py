import pandas as pd
from scipy.stats import norm, gamma

daf = pd.DataFrame({'CHR': [1] * 5,
                    'BP': range(1,6),
                    'REF': ['A'] * 5,
                    'ALT': ['C'] * 5,
                    'BETA': norm.rvs(size = 5),
                    'SE': gamma.rvs(a = 1, scale = 1, size = 5)})

daf.to_csv("resources/gwas/check-file-separator.tsv.gz", compression = 'gzip', index = False, sep = '\t')
daf.to_csv("resources/gwas/check-file-separator.csv.gz", compression = 'gzip', index = False, sep = ',')
daf.to_csv("resources/gwas/check-file-separator.semicolon.gz", compression = 'gzip', index = False, sep = ';')
daf.to_csv("resources/gwas/check-file-separator.space.gz", compression = 'gzip', index = False, sep = ' ')
