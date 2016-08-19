import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv


wd = '/Users/emanuel/Projects/projects/frezza_fh/'

samples = read_csv('%s/tables/kirp_samples.txt' % wd, sep='\t')

data = read_csv('/Users/emanuel/Projects/projects/tcga_rna_seq/data_preprocessed/KIRP_rsem_raw_counts.tab', sep='\t', index_col=0)
print '[INFO] data imported'

data_healthy = data[[c for c in data if c.split('-')[3][:-1] == '11']]

data = data.drop(set(data_healthy), axis=1)
data.columns = ['-'.join(c.split('-')[:3]) for c in data]

data_tumour = data[list(samples['code'])]
data_healthy.columns = ['-'.join(c.split('-')[:3]) for c in data_healthy]

data_healthy.to_csv('/Users/emanuel/Projects/projects/frezza_fh/tables/kirp_healthy.txt', sep='\t')
data_tumour.to_csv('/Users/emanuel/Projects/projects/frezza_fh/tables/kirp_tumour.txt', sep='\t')
