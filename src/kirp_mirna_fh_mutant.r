library(limma)
library(edgeR)

# -- FH mutants
fh_mutants <- c('TCGA.BQ.5879', 'TCGA.BQ.5893', 'TCGA.BQ.5894', 'TCGA.F9.A4JJ', 'TCGA.GL.7966')

# -- Import data
normal <- read.csv('/Users/emanuel/Projects/projects/frezza_fh/tables/kirp_mirna_healthy.txt', sep='\t', row.names=1)
tumour <- read.csv('/Users/emanuel/Projects/projects/frezza_fh/tables/kirp_mirna_tumour.txt', sep='\t', row.names=1)[fh_mutants]

# -- Design
design <- cbind(
  tumour = c(rep(1, dim(tumour)[2]), rep(0, dim(normal)[2])),
  normal = c(rep(0, dim(tumour)[2]), rep(1, dim(normal)[2]))
)

# -- Bind data-sets
dataset <- cbind(tumour, normal)

# -- Create a DGEList object using the edgeR package 
dataset <- DGEList(counts=dataset, genes=rownames(dataset))

# -- The limma-voom method assumes that rows with zero or very low counts have been removed.
dataset <- dataset[aveLogCPM(dataset) > 1,] # average CPM > 3

# -- Normalise with TMT
dataset <- calcNormFactors(dataset, method='TMM')

# -- Voom transformation
dataset <- voom(dataset)

# -- Limma differential expression
fit <- lmFit(dataset, design)
cont_matrix <- makeContrasts(tumourvsnormal=tumour-normal, levels=design)

fit_2 <- contrasts.fit(fit, cont_matrix)
fit_2 <- eBayes(fit_2)

# -- Retrieve results
res <- as.data.frame(topTable(fit_2, adjust.method='bonferroni', n=Inf))

# -- Export results
write.table(res, '/Users/emanuel/Projects/projects/frezza_fh/tables/kirp_mirna_limma_fh_mutant_vs_normal.txt', sep='\t', quote=F, row.names=F)
