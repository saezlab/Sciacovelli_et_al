import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
from pandas import read_csv, DataFrame, Series


def signif(value):
    if value < 0.01:
        return '*'
    elif value < 0.05:
        return '**'
    else:
        return '-'


def plot_volcanos_results(f, d, x_label, y_label, adj_p_value_label, title='', to_highlight=None, ids=None, type=None):
    sns.set(style='ticks')
    sns.set_context('paper')
    sns.despine(offset=10)

    # Descritise significance
    d['signif'] = [signif(v) for v in d[adj_p_value_label]]

    # Define pallete
    colour_pallete = [sns.color_palette('OrRd')[4], sns.color_palette('OrRd')[3], sns.color_palette('OrRd')[1]]
    sns.lmplot(x=x_label, y=y_label, data=d, hue='signif', fit_reg=False, palette=colour_pallete, legend=False)

    # Add FDR threshold lines
    plt.text(plt.xlim()[0]*.98, -np.log10(0.01), 'FDR 1%', ha='left', color=colour_pallete[0], alpha=0.65, fontsize=7)
    plt.axhline(-np.log10(0.01), c=colour_pallete[0], ls='--', lw=.5, alpha=.7)

    plt.text(plt.xlim()[0]*.98, -np.log10(0.05), 'FDR 5%', ha='left', color=colour_pallete[1], alpha=0.65, fontsize=6)
    plt.axhline(-np.log10(0.05), c=colour_pallete[1], ls='--', lw=.5, alpha=.7)

    # Add axis lines
    plt.axvline(0, c='#95a5a6', lw=.3, alpha=.15)
    plt.axhline(0, c='#95a5a6', lw=.3, alpha=.15)

    # Add axis labels and title
    plt.title(title, fontsize=10, fontname='Arial')
    plt.xlabel('fold change', fontsize=8, fontname='Arial')
    plt.ylabel('adj. p-value', fontsize=8, fontname='Arial')

    # Add text to highlighted genes
    if to_highlight is not None:
        for i, r in d.iterrows():
            if type is None:
                if r[ids] in to_highlight:
                    plt.text(r[x_label] * 1.01, r[y_label] * 1.01, r[ids], ha='left', alpha=0.75, fontsize=8)

            elif type == 'trans':
                # if r[adj_p_value_label] < 0.05:
                row_ids = set(gene.split('//')[1].strip() for gene in r['geneassignment'].split('///'))
                row_ids = row_ids.intersection(set(to_highlight))

                for row_id in row_ids:
                    plt.text(r[x_label] * 1.01, r[y_label] * 1.01, row_id, ha='left', alpha=0.75, fontsize=8)

    # Adjust axis lines thickness
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.xaxis.set_tick_params(width=0.5)
    ax.yaxis.set_tick_params(width=0.5)

    # Save plot
    fig = plt.gcf()
    fig.set_size_inches(5., 8.)
    fig.savefig(f)

    print '[INFO] Volcano generated: ' + f

data_dir = '/Users/emanuel/Projects/data/fh_cells/'
wd = '/Users/emanuel/Projects/projects/frezza_fh/'

# Import ensembl ids map
ensembl_human = read_csv('%s/tables/ensembl_id_map_human.txt' % wd, sep='\t', index_col=0).dropna().to_dict()['HGNC symbol']
ensembl_mouse = read_csv('%s/tables/ensembl_id_map_mouse.txt' % wd, sep='\t', index_col=0).dropna().to_dict()['MGI symbol']


# mi-RNAs
human_mirna = read_csv('%s/human_miRNAs/131122_human_microRNA.tsv' % data_dir, sep='\t', header=0).dropna(subset=['adj.P.Val', 'logFC'])
human_mirna['p.value.log10'] = -np.log10(human_mirna['adj.P.Val'])
human_mirna['ID'] = [id[4:] for id in human_mirna['ID']]
plot_volcanos_results(
    wd + 'report/marco_final/volcano_human_mirna.pdf',
    human_mirna,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - KO vs WT - miRNA',
    ['miR-200a', 'miR-200b', 'miR-200c', 'miR-141', 'miR-429', 'miR-486*', 'miR-1960', 'miR-1892', 'miR-34b-5p', 'miR-34c', 'miR-33', 'miR-136', 'miR-517b', 'miR-511', 'miR-196a'],
    'ID'
)
plt.close('all')


mouse_mirna = read_csv('%s/mouse_miRNAs/mouse_miRNAs.tsv' % data_dir, sep='\t', header=0).dropna(subset=['adj.P.Val', 'logFC'])
mouse_mirna['logFC'] *= -1
mouse_mirna['p.value.log10'] = -np.log10(mouse_mirna['adj.P.Val'])
mouse_mirna['ID'] = [i[4:] for i in mouse_mirna['ID']]
plot_volcanos_results(
    wd + 'report/marco_final/volcano_mouse_mirna.pdf',
    mouse_mirna,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Mouse - KO vs WT - miRNA',
    ['miR-200a', 'miR-200b', 'miR-200c', 'miR-141', 'miR-429', 'miR-486*', 'miR-1960', 'miR-1892', 'miR-34b-5p', 'miR-34c', 'miR-33'],
    'ID'
)
plt.close('all')


# KIRP mi-rna
kirp_mirna = read_csv('%s/tables/kirp_mirna_limma_tumour_vs_normal.txt' % wd, sep='\t')
kirp_mirna['p.value.log10'] = -np.log10(kirp_mirna['adj.P.Val'])
plot_volcanos_results(
    wd + 'report/marco_final/volcano_kirp_mirna.pdf',
    kirp_mirna,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'KIRP - Tumour vs Normal - miRNA',
    ['hsa-mir-200a', 'hsa-mir-200b', 'hsa-mir-200c', 'hsa-mir-141', 'hsa-mir-429'],
    'genes'
)
plt.close('all')


# KIRP mi-rna fh mutant
kirp_mirna = read_csv('%s/tables/kirp_mirna_limma_fh_mutant_vs_normal.txt' % wd, sep='\t')
kirp_mirna['p.value.log10'] = -np.log10(kirp_mirna['adj.P.Val'])
plot_volcanos_results(
    wd + 'report/marco_final/volcano_kirp_fh_mutant_mirna.pdf',
    kirp_mirna,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'KIRP - Tumour (fh mutant) vs Normal - miRNA',
    ['hsa-mir-200a', 'hsa-mir-200b', 'hsa-mir-200c', 'hsa-mir-141', 'hsa-mir-429'],
    'genes'
)
plt.close('all')


# Proteomics
human_proteomics = read_csv('%s/tables/protein_human_limma.tab' % wd, sep='\t', header=0).dropna(subset=['adj.P.Val', 'logFC'])
human_proteomics['ID'] = human_proteomics.index
plot_volcanos_results(
    wd + 'report/marco_final/volcano_human_proteomics.pdf',
    human_proteomics,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - KO vs WT - proteomics',
    ['VIM', 'FH'],
    'ID'
)
plt.close('all')

mouse_proteomics = read_csv('%s/tables/protein_mouse_limma.tab' % wd, sep='\t', header=0).dropna(subset=['adj.P.Val', 'logFC'])
mouse_proteomics['ID'] = mouse_proteomics.index
plot_volcanos_results(
    wd + 'report/marco_final/volcano_mouse_proteomics.pdf',
    mouse_proteomics,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Mouse - KO vs WT - proteomics',
    ['VIME_MOUSE', 'FH'],
    'ID'
)
plt.close('all')

# RNA-seq: Human UOK262
human_rnaseq = read_csv('%s/tables/human_rnaseq_UOK262_vs_UOK262pFH_limma.txt' % wd, sep='\t', index_col=0).dropna(subset=['adj.P.Val', 'logFC'])
human_rnaseq['ID'] = [ensembl_human[i] if i in ensembl_human else '' for i in human_rnaseq.index]
human_rnaseq['p.value.log10'] = -np.log10(human_rnaseq['adj.P.Val'])
plot_volcanos_results(
    '%s/report/marco_final/volcano_human_rnaseq_UOK262.pdf' % wd,
    human_rnaseq,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - KO vs WT - RNA-seq',
    ['VIM', 'FH', 'CDH1', 'CDH2'],
    'ID'
)
plt.close('all')

# RNA-seq: Mouse Clone1
mouse_rnaseq = read_csv('%s/tables/mouse_rnaseq_Clone1_vs_FH1_FL_FL_limma.txt' % wd, sep='\t', index_col=0).dropna(subset=['adj.P.Val', 'logFC'])
mouse_rnaseq['ID'] = [ensembl_mouse[i] if i in ensembl_mouse else '' for i in mouse_rnaseq.index]
mouse_rnaseq['p.value.log10'] = -np.log10(mouse_rnaseq['adj.P.Val'])
plot_volcanos_results(
    '%s/report/marco_final/volcano_mouse_rnaseq_clone1.pdf' % wd,
    mouse_rnaseq,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Mouse - Clone 1 vs FH fl/fl - RNA-seq',
    ['Vim', 'Fh1', 'Cdh1', 'Cdh2'],
    'ID'
)
plt.close('all')

# RNA-seq: Mouse Clone19
mouse_rnaseq = read_csv('%s/tables/mouse_rnaseq_Clone19_vs_FH1_FL_FL_limma.txt' % wd, sep='\t', index_col=0).dropna(subset=['adj.P.Val', 'logFC'])
mouse_rnaseq['ID'] = [ensembl_mouse[i] if i in ensembl_mouse else '' for i in mouse_rnaseq.index]
mouse_rnaseq['p.value.log10'] = -np.log10(mouse_rnaseq['adj.P.Val'])
plot_volcanos_results(
    '%s/report/marco_final/volcano_mouse_rnaseq_clone19.pdf' % wd,
    mouse_rnaseq,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Mouse - Clone 19 vs FH fl/fl - RNA-seq',
    ['Vim', 'Fh1', 'Cdh1', 'Cdh2'],
    'ID'
)
plt.close('all')

# RNA-seq: Mouse REC
mouse_rnaseq = read_csv('%s/tables/mouse_rnaseq_CLONE19_vs_REC_limma.txt' % wd, sep='\t', index_col=0).dropna(subset=['adj.P.Val', 'logFC'])
mouse_rnaseq['ID'] = [ensembl_mouse[i] if i in ensembl_mouse else '' for i in mouse_rnaseq.index]
mouse_rnaseq['p.value.log10'] = -np.log10(mouse_rnaseq['adj.P.Val'])
plot_volcanos_results(
    '%s/report/marco_final/volcano_mouse_rnaseq_Clone19_vs_REC.pdf' % wd,
    mouse_rnaseq,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Mouse - Clone 19 vs REC - RNA-seq',
    ['Vim', 'Fh1', 'Cdh1', 'Cdh2'],
    'ID'
)
plt.close('all')
