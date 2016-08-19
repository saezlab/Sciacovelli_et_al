import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv, melt

wd = '/Users/emanuel/Projects/projects/frezza_fh/'
data = '/Users/emanuel/Projects/data/fh_cells/'

mir = ['MIR200A', 'MIR200B', 'MIR429']

# -- Human methylation
h_data = read_csv('%s/human_methylation/UOK_pFh_vs_UOK_min_beta_annotated.txt' % data, sep='\t').dropna(subset=['UCSC_RefGene_Name'])
print 'Human 450k: ', h_data.shape

plot_df = h_data[[g in mir for g in h_data['UCSC_RefGene_Name']]]
plot_df['UOK 262 pFH'] = plot_df[[c for c in plot_df if c.startswith('UOK 262 pFH_min_')]].mean(1)
plot_df['UOK 262'] = plot_df[[c for c in plot_df if c.startswith('UOK 262_min')]].mean(1)
plot_df = plot_df[['UCSC_RefGene_Name', 'UOK 262 pFH', 'UOK 262']]
plot_df = melt(plot_df, id_vars=['UCSC_RefGene_Name'])

sns.set(style='ticks')
sns.boxplot('UCSC_RefGene_Name', 'value', 'variable', plot_df, order=mir, palette='Paired', sym='')
sns.stripplot('UCSC_RefGene_Name', 'value', 'variable', plot_df, order=mir, palette='Paired', jitter=True)
sns.despine(trim=True)
plt.xlabel('')
plt.ylabel('Methylation (b-value)')
plt.title('UOK 262')
plt.savefig('%s/report/marco_final/uok262_450k_methylation_mir.pdf' % wd, bbox_inches='tight')
plt.close('all')


# -- Mouse methylation
m_data = read_csv('%s/mouse_methylation/Fh1_fl_vs_Fh1_min_beta_annot.txt' % data, sep='\t').dropna(subset=['UCSC_REFGENE_NAME', 'UCSC_CPG_ISLANDS_NAME'])
# chr4:155451860-155452312
# chr4:155408917-155409306

m_data = m_data[[c.startswith('chr4') for c in m_data['UCSC_CPG_ISLANDS_NAME']]]
m_data[[int(c.split(':')[1].split('-')[0]) > (155408917 - 100) and int(c.split(':')[1].split('-')[1]) < (155452312 + 100) for c in m_data['UCSC_CPG_ISLANDS_NAME']]].T

m_data[[g in mir for g in m_data['UCSC_REFGENE_NAME']]]
