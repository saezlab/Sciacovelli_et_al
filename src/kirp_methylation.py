import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pandas import DataFrame, Series, read_csv

wd = '/Users/emanuel/Projects/projects/frezza_fh/'

mutations = read_csv('%s/tables/kirp_mutations.txt' % wd, sep='\t', index_col=0)
mutations['all'] = 1
mutations_dict = mutations.replace(0, np.NaN).unstack().reset_index().dropna().groupby('level_0')['code'].agg(lambda x: set(x)).to_dict()

data = read_csv('%s/tables/kirp_methylation_beta.tab' % wd, sep='\t', index_col=0)
data.columns = ['-'.join(c.split('-')[:3]) for c in data]

data = data[list(set(mutations.index).intersection(data))]

# -- Illumina 450k manifest
promoter_region = set(read_csv('/Users/emanuel/Projects/resources/illumina450k/cpg_sites_promoter_region.tab', sep='\t')['island'])

manifest = read_csv('/Users/emanuel/Projects/resources/illumina450k/HumanMethylation450_15017482_v1-2.csv', header=0, skiprows=7, index_col=0)
manifest = manifest.dropna(subset=['UCSC_RefGene_Name', 'UCSC_CpG_Islands_Name'])

print manifest


# --
mir = {'MIR200A', 'MIR200B', 'MIR200C'}
mir_manifest = manifest[[len(set(i.split(';')).intersection(mir)) > 0 for i in manifest['UCSC_RefGene_Name']]]
mir_manifest = mir_manifest[[i in promoter_region for i in mir_manifest['UCSC_CpG_Islands_Name']]]
mir_manifest = mir_manifest[['UCSC_RefGene_Name', 'UCSC_CpG_Islands_Name']]


# -- Plot
island_name = list(set(mir_manifest['UCSC_CpG_Islands_Name']))[0]

plot_df = DataFrame(data.ix[set(mir_manifest['UCSC_CpG_Islands_Name'].index)].mean(), columns=['beta'])
plot_df = [(k, m, v) for m in mutations.sum().sort(ascending=False, inplace=False).index for k, v in plot_df.ix[mutations_dict[m]].to_dict()['beta'].items()]
plot_df = DataFrame(plot_df, columns=['barcode', 'mutation', 'beta']).dropna()
plot_df['CpG island'] = island_name

sns.set(style='ticks')
g = sns.FacetGrid(plot_df, size=5, aspect=1)
g.map(sns.boxplot, 'beta', 'mutation', sym='', color='#95a5a6', orient='h')
g.map(sns.stripplot, 'beta', 'mutation', jitter=True, color='#95a5a6', orient='h')
g.map(plt.axhline, y=.5, lw=.3, ls='-', c='gray')
g.map(plt.axvline, x=.5, lw=.3, ls='--', c='gray')
sns.despine()
plt.title('KIRP')
plt.xlabel('promoter region of mir200a/b')
plt.ylabel('methylation (beta)')
plt.savefig('%s/report/marco_final/kirp_mir200_methylation_beta.pdf' % wd, bbox_inches='tight')
plt.close('all')
