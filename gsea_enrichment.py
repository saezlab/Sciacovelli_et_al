import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
from pymist.enrichment.gsea import gsea
from pandas import DataFrame, Series, read_csv, concat
from statsmodels.stats.multitest import multipletests
from matplotlib.gridspec import GridSpec


def gsea_barplot(dataset, filename):
    sns.set(style='ticks')
    plt.figure(figsize=(3, .3 * len(dataset)))

    colours, y_pos = [rgb2hex((r, g, b)) for r, g, b in sns.color_palette('Paired')[:2]], [x + 1.5 for x in range(len(dataset))]

    plt.barh(y_pos, -np.log10(dataset['pvalue']), lw=0, align='center', height=.5, color=colours[0], label='p-value')
    plt.barh(y_pos, -np.log10(dataset['fdr']), lw=0, align='center', height=.5, color=colours[1], label='FDR')
    plt.yticks(y_pos, dataset.index)

    plt.axvline(-np.log10(0.05), ls='--', lw=0.4, c='gray')
    plt.axvline(-np.log10(0.01), ls='--', lw=0.4, c='gray')

    plt.text(-np.log10(0.05) * 1.01, .5, '5%', ha='left', color='gray', fontsize=9)
    plt.text(-np.log10(0.01) * 1.01, .5, '1%', ha='left', color='gray', fontsize=9)

    sns.despine()
    plt.xlabel('-log10')
    plt.title('KEGG pathways enrichment')
    plt.legend(loc=0)
    plt.savefig('%s/report/marco_final/enrichments/%s.pdf' % (wd, filename), bbox_inches='tight')
    plt.close('all')
    print '[INFO] KEGG pathways enrichment plotted'


# -- Define variables
wd, n_permutations = '/Users/emanuel/Projects/projects/frezza_fh/', 10000


# -- Import gene signatures
# Import MSigDB kegg pathways
signatures = {}
with open('%s/tables/c2.cp.kegg.v5.0.symbols.gmt' % wd) as f:
    signatures = {' '.join(l.split('\t')[0][5:].lower().split('_')): set(l.strip().split('\t')[2:]) for l in f.readlines()}

# Import EMT signature
emt = read_csv('%s/tables/emt_human.tab' % wd, sep='\t')
signatures.update({s: set(emt.loc[emt['signature'] == s, 'gene']) for s in set(emt['signature'])})

# Import mouse/human orthologs
human2mouse = read_csv('%s/tables/human_mouse_orthologs.txt' % wd, sep='\t').dropna()
human2mouse = human2mouse.groupby('HGNC symbol').agg(lambda x: set(x)).to_dict()['MGI symbol']

# Convert human signatures gene symbols to mouse gene symbols
signatures_mouse = {s: {mg for hg in signatures[s] if hg in human2mouse for mg in human2mouse[hg]} for s in signatures}


# -- Import ensembl ids map
ensembl = read_csv('%s/tables/ensembl_id_map_human.txt' % wd, sep='\t', index_col=0).dropna().to_dict()['HGNC symbol']
ensembl_mouse = read_csv('%s/tables/ensembl_id_map_mouse.txt' % wd, sep='\t', index_col=0).dropna().to_dict()['MGI symbol']


# -- UOK262 enrichment
human_rnaseq = read_csv('%s/tables/human_rnaseq_UOK262_vs_UOK262pFH_limma.txt' % wd, sep='\t', index_col=0).dropna(subset=['adj.P.Val', 'logFC'])
human_rnaseq['ID'] = [ensembl[i] if i in ensembl else np.NaN for i in human_rnaseq.index]
human_rnaseq = human_rnaseq.dropna().groupby('ID')['t'].first().to_dict()

human_gsea = {s: gsea(human_rnaseq, signatures[s], n_permutations) for s in signatures}
human_gsea = DataFrame(human_gsea, index=['es', 'pvalue']).T
human_gsea['fdr'] = multipletests(human_gsea['pvalue'], method='fdr_bh')[1]
human_gsea = human_gsea.sort('fdr')

human_gsea.to_csv('%s/report/marco_final/enrichments/MSigDB_KEGG_UOK262.csv' % wd)

gsea_barplot(human_gsea[human_gsea['pvalue'] < 0.05].sort('fdr', ascending=False), 'MSigDB_KEGG_UOK262')

for signature in set(emt['signature']):
    plot = '%s/report/marco_final/enrichments/MSigDB_KEGG_UOK262_%s.pdf' % (wd, signature)
    title = 'UOK262 - %s' % signature
    gsea(human_rnaseq, set(emt.loc[emt['signature'] == signature, 'gene']), n_permutations, plot_name=plot, plot_title=title)

print '[INFO] GSEA done'


# -- SDH mutant cells
fc_ko = read_csv('%s/tables/FC_annot_eBayes_RMAnorm_SC8_flox_KO5_KO7.txt' % wd, sep='\t', index_col=0)[['Gene Symbol', 'KO5vsFlox', 'KO7vsFlox']].dropna()
fc_ko = fc_ko.groupby('Gene Symbol').median()
fc_ko = fc_ko.to_dict()

for condition in set(fc_ko):
    mouse_gsea = {s: gsea(fc_ko[condition], signatures_mouse[s], n_permutations) for s in signatures_mouse}
    mouse_gsea = DataFrame(mouse_gsea, index=['es', 'pvalue']).T
    mouse_gsea = mouse_gsea.dropna()
    mouse_gsea['fdr'] = multipletests(mouse_gsea['pvalue'], method='fdr_bh')[1]
    mouse_gsea = mouse_gsea.sort('fdr')

    mouse_gsea.to_csv('%s/report/marco_final/enrichments/MSigDB_KEGG_%s.csv' % (wd, condition))

    gsea_barplot(mouse_gsea[mouse_gsea['pvalue'] < 0.05].sort('fdr', ascending=False), 'MSigDB_KEGG_%s.pdf' % condition)

    for signature in set(emt['signature']):
        plot = '%s/report/marco_final/enrichments/MSigDB_KEGG_%s_%s.pdf' % (wd, condition, signature)
        title = '%s - %s' % (condition, signature)
        gsea(fc_ko[condition], signatures_mouse[signature], n_permutations,  plot_name=plot, plot_title=title)

print '[INFO] GSEA done'


# --- Mouse: Clone1/19 vs FH1_FL_FL
for condition in ['Clone1', 'Clone19']:
    # Import RNA-seq data
    mouse_rnaseq = read_csv('%s/tables/mouse_rnaseq_%s_vs_FH1_FL_FL_limma.txt' % (wd, condition), sep='\t', index_col=0).dropna(subset=['adj.P.Val', 'logFC'])
    mouse_rnaseq['ID'] = [ensembl_mouse[i] if i in ensembl_mouse else np.NaN for i in mouse_rnaseq.index]
    mouse_rnaseq = mouse_rnaseq.dropna().groupby('ID')['t'].first().to_dict()

    mouse_gsea = {s: gsea(mouse_rnaseq, signatures_mouse[s], n_permutations) for s in signatures_mouse}
    mouse_gsea = DataFrame(mouse_gsea, index=['es', 'pvalue']).T
    mouse_gsea['fdr'] = multipletests(mouse_gsea['pvalue'], method='fdr_bh')[1]
    mouse_gsea = mouse_gsea.sort('fdr')
    mouse_gsea.to_csv('%s/report/marco_final/enrichments/MSigDB_KEGG_%s_vs_FH1_FL_FL.csv' % (wd, condition))

    gsea_barplot(mouse_gsea[mouse_gsea['pvalue'] < 0.05].sort('fdr', ascending=False), 'MSigDB_KEGG_%s_vs_FH1_FL_FL.pdf' % condition)

    for signature in set(emt['signature']):
        plot = '%s/report/marco_final/enrichments/MSigDB_KEGG_%s_vs_FH1_FL_FL_%s.pdf' % (wd, condition, signature)
        title = '%s vs FH1_FL_FL - %s' % (condition, signature)
        gsea(mouse_rnaseq, signatures_mouse[signature], n_permutations,  plot_name=plot, plot_title=title)

print '[INFO] GSEA done'


# ---- Mouse EMT: Clone19 vs REC
mouse_rnaseq = read_csv('%s/tables/mouse_rnaseq_CLONE19_vs_REC_limma.txt' % wd, sep='\t', index_col=0).dropna(subset=['adj.P.Val', 'logFC'])
mouse_rnaseq['ID'] = [ensembl_mouse[i] if i in ensembl_mouse else np.NaN for i in mouse_rnaseq.index]
mouse_rnaseq = mouse_rnaseq.dropna().groupby('ID')['t'].first().to_dict()

mouse_gsea = {s: gsea(mouse_rnaseq, signatures_mouse[s], n_permutations) for s in signatures_mouse}
mouse_gsea = DataFrame(mouse_gsea, index=['es', 'pvalue']).T
mouse_gsea['fdr'] = multipletests(mouse_gsea['pvalue'], method='fdr_bh')[1]
mouse_gsea = mouse_gsea.sort('fdr')
mouse_gsea.to_csv('%s/report/marco_final/enrichments/MSigDB_KEGG_CLONE19_vs_REC.csv' % wd)

gsea_barplot(mouse_gsea[mouse_gsea['pvalue'] < 0.05].sort('fdr', ascending=False), 'MSigDB_KEGG_CLONE19_vs_REC.pdf')

for signature in set(emt['signature']):
    plot = '%s/report/marco_final/enrichments/MSigDB_KEGG_CLONE19_vs_REC_%s.pdf' % (wd, signature)
    title = '%s - %s' % ('CLONE19 vs REC', signature)
    gsea(mouse_rnaseq, signatures_mouse[signature], n_permutations,  plot_name=plot, plot_title=title)

print '[INFO] GSEA done'


# -- KIRP enrichment
kirp = read_csv('%s/tables/kirp_limma_tumour_vs_normal.txt' % wd, sep='\t', index_col=0)['t'].to_dict()

kirp_gsea = {s: gsea(kirp, signatures[s], n_permutations) for s in signatures}
kirp_gsea = DataFrame(kirp_gsea, index=['es', 'pvalue']).T
kirp_gsea['fdr'] = multipletests(kirp_gsea['pvalue'], method='fdr_bh')[1]
kirp_gsea = kirp_gsea.sort('fdr')
kirp_gsea.to_csv('%s/report/marco_final/enrichments/MSigDB_KEGG_KIRP_Tumour_vs_Normal.csv' % wd)
print '[INFO] GSEA done'

gsea_barplot(kirp_gsea[kirp_gsea['pvalue'] < 0.05].sort('fdr', ascending=False), 'MSigDB_KEGG_KIRP_Tumour_vs_Normal.pdf')

for signature in set(emt['signature']):
    plot = '%s/report/marco_final/enrichments/MSigDB_KEGG_KIRP_Tumour_vs_Normal_%s.pdf' % (wd, signature)
    title = '%s - %s' % ('KIRP (Tumour vs Normal)', signature)
    gsea(kirp, signatures[signature], n_permutations,  plot_name=plot, plot_title=title)

print '[INFO] GSEA done'


# -- KIRP enrichment FH mutant
kirp = read_csv('%s/tables/kirp_limma_fhmutant_vs_normal.txt' % wd, sep='\t', index_col=0)['t'].to_dict()

kirp_gsea = {s: gsea(kirp, signatures[s], n_permutations) for s in signatures}
kirp_gsea = DataFrame(kirp_gsea, index=['es', 'pvalue']).T
kirp_gsea['fdr'] = multipletests(kirp_gsea['pvalue'], method='fdr_bh')[1]
kirp_gsea = kirp_gsea.sort('fdr')
kirp_gsea.to_csv('%s/report/marco_final/enrichments/MSigDB_KEGG_KIRP_Tumour_FH_mutant_vs_Normal.csv' % wd)
print '[INFO] GSEA done'

gsea_barplot(kirp_gsea[kirp_gsea['pvalue'] < 0.05].sort('fdr', ascending=False), 'MSigDB_KEGG_KIRP_Tumour_FH_mutant_vs_Normal.pdf')

for signature in set(emt['signature']):
    plot = '%s/report/marco_final/enrichments/MSigDB_KEGG_KIRP_Tumour_FH_mutant_vs_Normal%s.pdf' % (wd, signature)
    title = '%s - %s' % ('KIRP FH mutant (Tumour vs Normal)', signature)
    gsea(kirp, signatures[signature], n_permutations,  plot_name=plot, plot_title=title)

print '[INFO] GSEA done'
