library(sva)
library(limma)
library(pheatmap)
library(MatrixGenerics)
library(RColorBrewer)
library(ggrepel)
library(PCAtools)

# Detect and model batch effects using SVA and limma
setwd(<working_dir>)
load("expression_data/GSE93776_affy_norm_outliers_removed.RData")

eset <- Affy_norm
colnames(eset)
colnames(exprs(eset)) <- rownames(pData(eset))
colnames(exprs(eset))

# Examine metadata
head(pData(eset))

# Use SVA to identify latent sources of variance
pheno <- pData(eset)
head(pheno)

# Create null model matrix with adjustment variables only
mod0 <- model.matrix(~ as.factor(Cell_type) + 
                       as.factor(Gender), 
                     data = pheno)

# Create full model matrix that includes variable of interest and
# adjustment variables
mod <- model.matrix(~ as.factor(Trait) + as.factor(Cell_type) + 
                       as.factor(Gender), 
                     data = pheno)

# Detect latent batch effects
edata <- exprs(eset)
svobj <- sva(edata, mod, mod0)
head(svobj$sv)

# Use limma to correct for hidden batch effects detected with SVA
# Specify design matrix for the variables we want to preserve
mod <- model.matrix(~ as.factor(Trait) + as.factor(Cell_type) +
                      as.factor(Gender), 
                    data = pheno)

# Create a matrix of batch corrected values
edata_cor <- removeBatchEffect(edata, covariates = svobj$sv, 
                               design = mod)
head(edata)[,1:3]
head(edata_cor)[,1:3]

# Cluster data before and after correction
# 1. Before correction
# Select top 500 genes with highest Mean Absolute Deviations (MAD)
mads_sorted <- sort(rowMads(edata), decreasing = T)
mads_sorted <- mads_sorted[1:500]
head(mads_sorted)

# Extract the matrix with top 500 most variable values
top500 <- edata[which(rownames(edata) %in% names(mads_sorted)),]

# Create column annotation for heatmap
rownames(sample_data) <- colnames(edata)
annot <- data.frame(row.names = rownames(pData(eset)),
                    Trait = pData(eset)$Trait,
                    Gender = pData(eset)$Gender,
                    Age_group = pData(eset)$Cell_type)

# Cluster genes and samples, plot heatmap
pdf("GSE93776_heatmap_before_correction.pdf", width = 7, height = 10)
pheatmap(top500, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annot, 
         show_rownames = F)
dev.off()

# ----------------------------Batch removal ---------------------------------- #
# After batch removal
mads_sorted <- sort(rowMads(edata_cor), decreasing = T)
mads_sorted <- mads_sorted[1:500]
head(mads_sorted)

# Extract the matrix with top 500 most variable values
top500 <- edata_cor[which(rownames(edata_cor) %in% names(mads_sorted)),]

# Create column annotation for heatmap
rownames(sample_data) <- colnames(edata)
annot <- data.frame(row.names = rownames(pData(eset)),
                    Trait = pData(eset)$Trait,
                    Gender = pData(eset)$Gender,
                    Cell_type = pData(eset)$Cell_type)

# Cluster genes and samples, plot heatmap
# Note that the samples are now clearly clustered according to the Trait
pdf("GSE93776_heatmap_after_correction.pdf", width = 7, height = 8)
pheatmap(top500, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annot, 
         show_rownames = F)
dev.off()

# -------------------------- PCA ---------------------------------------- #
# Plot PCA after batch correction
p <- pca(edata_cor, metadata = annot, removeVar = 0.9)

# Create biplot and draw stat ellipses at 95% CI around groups
biplot(p,
       colby = 'Cell_type', 
       #colkey = c('Control' = 'forestgreen', 'RA' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', 
       legendLabSize = 8, 
       legendIconSize = 4)
ggsave("93776_PCA_biplot_after_correct.pdf", device = "pdf", units = "in", 
       width = 6, height = 6)

# Create pairs plot
pairsplot(p,
          components = getComponents(p, c(1:6)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Trait',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
ggsave("GSE93776_pairs_plot_after_correct.pdf", device = "pdf", units = "in", 
       width = 6, height = 6)



# Plot loadings, fancy
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)
ggsave("GSE93776_loadings_plot_corrected.pdf", device = "pdf", units = "in", 
       width = 6, height = 7)

# Create eigencorrelation plot
pdf("GSE93776_eigencor_plot_correct.pdf", width = 7, height = 7) 
eigencorplot(p,
             metavars = c('Trait', 'Cell_type', 'Gender'))
dev.off()

# Identify optimal number of PCs
#horn <- parallelPCA(expr) # Horn's method
#horn$n # 4

elbow <- findElbowPoint(p$variance)
elbow # 3

# Create screeplot
screeplot(p, components = getComponents(p, 1:10),
          vline = c(elbow)) + 
  #geom_label(aes(x = horn$n + 1, y = 50,
  #               label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))
ggsave("GSE93776_screeplot.pdf", device = "pdf", units = "in", 
       width = 7, height = 7)

# -------------------- Detect DEGs after correction ----------------------- #
# Create expression set object with VSN normalized data
# Note, that we should not use SVA corrected data in differential expression
# analysis. Instead, we should supply latent variables detected with SVA
# as covariates in model matrix
surrogate <- svobj$sv
surrogate <- surrogate[,1:7]
colnames(surrogate) <- c("SV1", "SV2", "SV3", "SV4", "SV5",
                         "SV6", "SV7")
# Add surrogate variables to pheno data
pheno <- cbind(pData(eset), surrogate) 
head(pheno)
# Create model matrix with adjustment and surrogate variables
design <- model.matrix(~ Trait + Cell_type + Gender + SV1 +
                         SV2 + SV3 + SV4 + SV5 + SV6 + SV7, 
                       data = pheno)
head(design)

# Fit the model
fit <- lmFit(eset, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "TraitRA", number = Inf, 
                adjust.method = "BH") 
head(res)
write.csv(res, file = "GSE93776_RA_vs_Controls_corrected.csv")

# How many genes are significant (adj.p.value < 0.05 and absolute 
# fold change > 1.5)
sum(res$adj.P.Val < 0.05) # 5545 probes
# Save the results for significant probes
write.csv(res[res$adj.P.Val < 0.05,], 
          file = "GSE93776_RA_vs_Controls_sig_only_corrected.csv")

# Plot a heatmap of differentially expressed probes
sig_probes <- res[which(res$adj.P.Val < 0.05 &
                          abs(res$logFC) > 0.59),]
sig_probes <- rownames(sig_probes)
expr <- edata_cor
sig_exp <- expr[rownames(expr) %in% sig_probes,]

annot <- data.frame(Trait = pData(eset)$Trait,
                    Gender = pData(eset)$Gender,
                    Cell_type = pData(eset)$Cell_type)
rownames(annot) <- rownames(pData(eset))
head(annot)

# Clustering and heatmap for significant probes (adj. p-values < 0.05)
pdf("GSE93776_heatmap_significant_corrected.pdf", width = 7, height = 6)
pheatmap(sig_exp,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annot,
         show_rownames = F)
dev.off()

# Create MA plot showing the relationship between average expression
# and log fold change
res <- res[with(res, order(logFC)),]
res$threshold <- as.factor(res$adj.P.Val < 0.05)
ggplot(data=res, aes(x=log2(AveExpr), y=logFC, colour=threshold)) +
  geom_point(alpha=0.4, size=1.8) +
  geom_hline(aes(yintercept = 0), colour = "blue", size = 1.2) +
  ylim(c(min(res$logFC), max(res$logFC))) +
  xlab("Mean expression") +
  ylab("Log2 Fold Change") +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 12)) +
  theme(axis.title.y = element_text(face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 12)) +
  scale_colour_discrete(name = "p.adjusted < 0.05") +
  theme(legend.title = element_text(face = "bold", size = 15)) +
  theme(legend.text = element_text(size = 14))
ggsave("GSE93776_MAplot_corrected.pdf", device = "pdf", width = 7, height = 5)

# Create volcano plot
res <- res[with(res, order(logFC)),]
#res$threshold <- as.factor(abs(res$logFC) > 0.59 & res$adj.P.Val < 0.05)
res$threshold <- as.factor(res$adj.P.Val < 0.05)
ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) + xlim(c(min(res$logFC),
                                            c(max(res$logFC)))) +
  ylim(c(min(-log10(res$adj.P.Val)), max(-log10(res$adj.P.Val)))) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(axis.text=element_text(size=12, face="bold")) +
  theme(axis.title=element_text(size=14)) +
  theme(legend.title=element_text(size=14)) +
  theme(legend.text=element_text(size=12))
ggsave("GSE93776_volcano_plot_corrected.pdf", device = "pdf",
       width = 5, height = 7)
