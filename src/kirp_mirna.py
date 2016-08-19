import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv


wd = '/Users/emanuel/Projects/projects/frezza_fh/'
samples = list(read_csv('%s/tables/kirp_samples.txt' % wd, sep='\t')['code'])

# fh_mutants = ['TCGA-BQ-5879', 'TCGA-BQ-5893', 'TCGA-BQ-5894', 'TCGA-F9-A4JJ', 'TCGA-GL-7966']
# fh_mutants = ['TCGA-BQ-5879-01A-11R-1591-13', 'TCGA-BQ-5893-01A-11R-1591-13', 'TCGA-BQ-5894-01A-11R-1591-13', 'TCGA-F9-A4JJ-01A-11R-A25B-13', 'TCGA-GL-7966-01A-11R-2203-13']
cimp = ['TCGA-A4-7915-01A-11R-2203-13', 'TCGA-BQ-5879-01A-11R-1591-13', 'TCGA-BQ-5893-01A-11R-1591-13', 'TCGA-BQ-5894-01A-11R-1591-13', 'TCGA-F9-A4JJ-01A-11R-A25B-13', 'TCGA-G7-6793-01A-11R-1964-13', 'TCGA-GL-7966-01A-11R-2203-13', 'TCGA-P4-A5E8-01A-11R-A28L-13', 'TCGA-P4-A5EA-01A-11R-A28L-13']

# -- Subset
data = read_csv('%s/tables/kirp_mirnaseq.txt' % wd, sep='\t', index_col=0)
data = DataFrame(data.loc[:, data.ix['miRNA_ID'] == 'read_count'].drop('miRNA_ID'), dtype=float)

data_healthy = data[[c for c in data if (c.split('-')[3][:-1] == '11') and ('-'.join(c.split('-')[:3]) in samples)]]

data = data.drop(set(data_healthy), axis=1)

data_tumour = data[[c for c in data if '-'.join(c.split('-')[:3]) in samples]]

data_healthy.to_csv('/Users/emanuel/Projects/projects/frezza_fh/tables/kirp_mirna_healthy.txt', sep='\t')
data_tumour.to_csv('/Users/emanuel/Projects/projects/frezza_fh/tables/kirp_mirna_tumour.txt', sep='\t')


res = read_csv('%s/tables/kirp_mirna_limma_fh_mutant_vs_normal.txt' % wd, sep='\t', index_col=0)
res.sort('logFC')

mirna = 'hsa-mir-200a'


# -- RMSE
data = read_csv('%s/tables/kirp_mirnaseq.txt' % wd, sep='\t', index_col=0)
data = data.loc[:, data.ix['miRNA_ID'] == 'reads_per_million_miRNA_mapped'].drop('miRNA_ID')
data.columns = [c.split('.')[0] for c in data]

# mi-rnas boxplots
mirnas = ['hsa-mir-200a', 'hsa-mir-200b', 'hsa-mir-200c', 'hsa-mir-141', 'hsa-mir-429']

plotdf = DataFrame(data.ix[mirnas], dtype=float).unstack().reset_index()
plotdf.columns = ['id', 'var', 'val']

plotdf['type'] = 'tumour'
plotdf.loc[[c in cimp for c in plotdf['id']], 'type'] = 'CIMP'
plotdf.loc[[c.split('-')[3][:-1] == '11' for c in plotdf['id']], 'type'] = 'normal'
plotdf = plotdf[['-'.join(c.split('-')[:3]) in samples for c in plotdf['id']]]

sns.set(style='white')
g = sns.FacetGrid(plotdf, col='var', sharey=False)
g.map(sns.boxplot, 'type', 'val', palette='Set1', order=['normal', 'tumour', 'CIMP'], sym='')
g.map(sns.stripplot, 'type', 'val', palette='Set1', order=['normal', 'tumour', 'CIMP'], jitter=True)
g.set_axis_labels('', 'reads per million')
plt.savefig('%s/report/marco_final/KIRP_mirnas_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')

#
mirnas = ['hsa-mir-200a', 'hsa-mir-200b', 'hsa-mir-200c', 'hsa-mir-141', 'hsa-mir-429']

plotdf = DataFrame(data.ix[mirnas], dtype=float).unstack().reset_index()
plotdf.columns = ['id', 'var', 'val']

plotdf['type'] = 'tumour'

plotdf = plotdf[['-'.join(c.split('-')[:3]) in fh_mutants for c in plotdf['id']]]
plotdf.loc[[c.split('-')[3][:-1] == '11' for c in plotdf['id']], 'type'] = 'normal'
