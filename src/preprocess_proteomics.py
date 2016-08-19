import numpy as np
from pandas import read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename

# Configure vars
data_dir = '/Users/emanuel/Projects/data/fh_cells/'
wd_dir = '/Users/emanuel/Projects/projects/frezza_fh/'

for organism in ['human', 'mouse']:
    os = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}[organism]

    # Import Uniprot gene name
    uniprot2genename = read_uniprot_genename(os=os)

    # Define files paths
    tp_file = {
        'human': 'human_proteomics/b1368p100_protein_human.tab',
        'mouse': 'mouse_proteomics/b1368p100_protein_mouse.tab'
    }[organism]

    tp_file_processed = {
        'human': 'tables/protein_human_processed.tab',
        'mouse': 'tables/protein_mouse_processed.tab'
    }[organism]

    # Import samplesheet
    ss = read_csv(data_dir + 'fh_samplesheet.tab', sep='\t', index_col=0)
    ss = ss.loc[ss['organism'] == organism]

    ss_tp = list(ss[ss['type'] == 'tp'].index)
    ss_tp_ko = list(ss[(ss['type'] == 'tp') & (ss['condition'] == 'fh_ko')].index)
    ss_tp_wt = list(ss[(ss['type'] == 'tp') & (ss['condition'] == 'fh_wt')].index)

    # ---- Import and preprocess phospho and proteomics data
    tp_all = read_csv(data_dir + tp_file, sep='\t').dropna(subset=['uniprot'])

    # Drop NaN on phospho sites and peptides
    tp_all = tp_all.dropna(subset=['peptide'])

    # Assemble reduced data-set
    tp_columns = ['uniprot'] if organism != 'mouse' else ['uniprot', 'acc_no']
    tp_columns.extend(ss_tp)
    tp = tp_all[tp_columns]

    # ---- Replace zeros with NaN
    tp = tp.replace(0.0, np.NaN)

    # ---- Remove peptides with 1 or less measurements per condition
    tp = tp.loc[(np.isnan(tp.ix[:, ss_tp_ko]).sum(axis=1) < 2) & (np.isnan(tp.ix[:, ss_tp_wt]).sum(axis=1) < 2), ]

    # ---- Remove ambigous peptides
    tp = tp[[len(x.split('; ')) == 2 and x != '; ' for x in tp['uniprot']]]
    tp['uniprot'] = [x.split('; ')[0] for x in tp['uniprot' if organism != 'mouse' else 'acc_no']]

    # ---- Average peptides to protein level
    tp = tp.groupby(by='uniprot').median()

    # ---- Log 2 transform
    tp[ss_tp] = np.log2(tp[ss_tp])

    # ---- Scale replicates
    tp[ss_tp] = (tp[ss_tp] - tp[ss_tp].mean()) / tp[ss_tp].std()

    # ---- Export
    tp.columns = [ss.ix[c, 'condition'] for c in tp.columns]
    tp.index = [uniprot2genename[c][0] if c in uniprot2genename else c for c in tp.index]

    tp.to_csv(wd_dir + tp_file_processed, sep='\t')

    # ---- Verbose
    print '[INFO] Preprocess done!'