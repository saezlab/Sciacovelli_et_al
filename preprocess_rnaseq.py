import os
import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from pandas import DataFrame, Series, read_csv, concat


def pearson(x, y):
    mask = np.bitwise_and(np.isfinite(x), np.isfinite(y))
    cor, pvalue = pearsonr(x[mask], y[mask]) if np.sum(mask) > 1 else (np.NaN, np.NaN)
    return cor, pvalue

data_dirs = [
    '/Users/emanuel/Projects/data/fh_cells/human_rnaseq/human/',
    '/Users/emanuel/Projects/data/fh_cells/mouse_rnaseq/mouse/'
]

wd = '/Users/emanuel/Projects/projects/frezza_fh/'

# Import RNA-seq data into a single matrix
human_gexp, mouse_gexp = [concat([read_csv(d + f, sep='\t', header=None, names=['id', f.split('.')[0]], index_col=0) for f in os.listdir(d)], axis=1) for d in data_dirs]
print '[INFO] Human/Mouse RNA-seq imported'

# Discard reads with techincal problems
human_gexp = human_gexp.ix[[i for i in human_gexp.index if not i.startswith('__')]]
mouse_gexp = mouse_gexp.ix[[i for i in mouse_gexp.index if not i.startswith('__')]]

# Export matricies
human_gexp.to_csv('%s/tables/human_rnaseq.txt' % wd, sep='\t')
mouse_gexp.to_csv('%s/tables/mouse_rnaseq.txt' % wd, sep='\t')
print '[INFO] Raw count tables exported'
