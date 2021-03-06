library(limma)

# Files paths
organism = 2 # 1: human, 2: mouse

tp_human_file = c(
  '/Users/emanuel/Projects/projects/frezza_fh/tables/protein_human_processed.tab', 
  '/Users/emanuel/Projects/projects/frezza_fh/tables/protein_mouse_processed.tab'
)[organism]

tp_result_file = c(
  '/Users/emanuel/Projects/projects/frezza_fh/tables/protein_human_limma.tab',
  '/Users/emanuel/Projects/projects/frezza_fh/tables/protein_mouse_limma.tab'
)[organism]

# Import data-sets
tp <- read.table(tp_human_file, sep='\t', row.names=1, header=T, stringsAsFactors=F, check.names=F)

# ---- Proteomics
# Run differential analysis
tp_design <- cbind(KO=rep(c(1,0), each=3), WT=rep(c(0,1), each=3))
tp_fit <- lmFit(tp, tp_design)
tp_cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=tp_design)
tp_fit_2 <- contrasts.fit(tp_fit, tp_cont_matrix)
tp_fit_2 <- eBayes(tp_fit_2)

tp_result <- as.data.frame(topTable(tp_fit_2, adjust.method='fdr', n=Inf))
tp_result$p.value.log10 <- -log10(tp_result$adj.P.Val)

# Store reults
write.table(tp_result, tp_result_file, sep='\t', quote=F)
