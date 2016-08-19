library(limma)
library(edgeR)

# ---- Human
# Import raw counts matrix - rows: genes, cols: samples
gexp <- read.csv('/Users/emanuel/Projects/projects/frezza_fh/tables/human_rnaseq.txt', sep='\t', row.names=1)
gexp <- gexp[, unlist(lapply(colnames(gexp), function(x) { grepl('^UOK262', x) & !grepl('SHRNA', x) }))]

# Create a DGEList object using the edgeR package 
gexp <- DGEList(counts=gexp,genes=rownames(gexp))

# The limma-voom method assumes that rows with zero or very low counts have been removed.
gexp <- gexp[aveLogCPM(gexp) > 1,] # average CPM > 3

# Normalise with TMT
gexp <- calcNormFactors(gexp, method='TMM')

# Voom transformation
gexp <- voom(gexp)

# Create UOK262 vs UOK262pFH design matrix
design <- cbind(
  UOK262pFH=as.integer(unlist(lapply(colnames(gexp), function(x) { grepl('^UOK262_PFH', x)}))),
  UOK262=as.integer(!unlist(lapply(colnames(gexp), function(x) { grepl('^UOK262_PFH', x)})))
)

# Limma differential expression
fit <- lmFit(gexp, design)
cont_matrix <- makeContrasts(UOK262vsUOK262pFH=UOK262-UOK262pFH, levels=design)

fit_2 <- contrasts.fit(fit, cont_matrix)
fit_2 <- eBayes(fit_2)

# Retrieve results
res <- as.data.frame(topTable(fit_2, adjust.method='bonferroni', n=Inf))

# Export results
write.table(res, '/Users/emanuel/Projects/projects/frezza_fh/tables/human_rnaseq_UOK262_vs_UOK262pFH_limma.txt', sep='\t', quote=F, row.names=F)


# ---- Mouse: Clone 1
# Import raw counts matrix - rows: genes, cols: samples
gexp <- read.csv('/Users/emanuel/Projects/projects/frezza_fh/tables/mouse_rnaseq.txt', sep='\t', row.names=1)
gexp <- gexp[, unlist(lapply(colnames(gexp), function(x) { grepl('^FH1_FL_FL_', x) || grepl('^CLONE1_', x) }))]

# Create a DGEList object using the edgeR package 
gexp <- DGEList(counts=gexp,genes=rownames(gexp))

# The limma-voom method assumes that rows with zero or very low counts have been removed.
gexp <- gexp[aveLogCPM(gexp) > 1,] # average CPM > 3

# Normalise with TMT
gexp <- calcNormFactors(gexp, method='TMM')

# Voom transformation
gexp <- voom(gexp)

# Create UOK262 vs UOK262pFH design matrix
design <- cbind(
  KO=as.integer(unlist(lapply(colnames(gexp), function(x) { grepl('^CLONE1_', x)}))),
  WT=as.integer(unlist(lapply(colnames(gexp), function(x) { grepl('^FH1_FL_FL_', x)})))
)

# Limma differential expression
fit <- lmFit(gexp, design)
cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=design)

fit_2 <- contrasts.fit(fit, cont_matrix)
fit_2 <- eBayes(fit_2)

# Retrieve results
res <- as.data.frame(topTable(fit_2, adjust.method='bonferroni', n=Inf))

# Export results
write.table(res, '/Users/emanuel/Projects/projects/frezza_fh/tables/mouse_rnaseq_Clone1_vs_FH1_FL_FL_limma.txt', sep='\t', quote=F, row.names=F)

# ---- Mouse: Clone 19
# Import raw counts matrix - rows: genes, cols: samples
gexp <- read.csv('/Users/emanuel/Projects/projects/frezza_fh/tables/mouse_rnaseq.txt', sep='\t', row.names=1)
gexp <- gexp[, unlist(lapply(colnames(gexp), function(x) { grepl('^FH1_FL_FL_', x) || grepl('^CLONE19_', x) }))]

# Create a DGEList object using the edgeR package 
gexp <- DGEList(counts=gexp,genes=rownames(gexp))

# The limma-voom method assumes that rows with zero or very low counts have been removed.
gexp <- gexp[aveLogCPM(gexp) > 1,] # average CPM > 3

# Normalise with TMT
gexp <- calcNormFactors(gexp, method='TMM')

# Voom transformation
gexp <- voom(gexp)

# Create UOK262 vs UOK262pFH design matrix
design <- cbind(
  KO=as.integer(unlist(lapply(colnames(gexp), function(x) { grepl('^CLONE19_', x)}))),
  WT=as.integer(unlist(lapply(colnames(gexp), function(x) { grepl('^FH1_FL_FL_', x)})))
)

# Limma differential expression
fit <- lmFit(gexp, design)
cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=design)

fit_2 <- contrasts.fit(fit, cont_matrix)
fit_2 <- eBayes(fit_2)

# Retrieve results
res <- as.data.frame(topTable(fit_2, adjust.method='bonferroni', n=Inf))

# Export results
write.table(res, '/Users/emanuel/Projects/projects/frezza_fh/tables/mouse_rnaseq_Clone19_vs_FH1_FL_FL_limma.txt', sep='\t', quote=F, row.names=F)


# ---- Mouse: REC
# Import raw counts matrix - rows: genes, cols: samples
gexp <- read.csv('/Users/emanuel/Projects/projects/frezza_fh/tables/mouse_rnaseq.txt', sep='\t', row.names=1)
gexp <- gexp[, unlist(lapply(colnames(gexp), function(x) { grepl('^REC_', x) || grepl('^CLONE19_', x) }))]

# Create a DGEList object using the edgeR package 
gexp <- DGEList(counts=gexp,genes=rownames(gexp))

# The limma-voom method assumes that rows with zero or very low counts have been removed.
gexp <- gexp[aveLogCPM(gexp) > 1,] # average CPM > 3

# Normalise with TMT
gexp <- calcNormFactors(gexp, method='TMM')

# Voom transformation
gexp <- voom(gexp)

# Create UOK262 vs UOK262pFH design matrix
design <- cbind(
  REC=as.integer(unlist(lapply(colnames(gexp), function(x) { grepl('^REC_', x)}))),
  CLONE19=as.integer(unlist(lapply(colnames(gexp), function(x) { grepl('^CLONE19_', x)})))
)

# Limma differential expression
fit <- lmFit(gexp, design)
cont_matrix <- makeContrasts(CLONE19vsREC=CLONE19-REC, levels=design)

fit_2 <- contrasts.fit(fit, cont_matrix)
fit_2 <- eBayes(fit_2)

# Retrieve results
res <- as.data.frame(topTable(fit_2, adjust.method='bonferroni', n=Inf))

# Export results
write.table(res, '/Users/emanuel/Projects/projects/frezza_fh/tables/mouse_rnaseq_CLONE19_vs_REC_limma.txt', sep='\t', quote=F, row.names=F)
