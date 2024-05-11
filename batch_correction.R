library(sva)
library(limma)
library(pheatmap)
library(MatrixGenerics)
library(RColorBrewer)
library(ggrepel)
library(PCAtools)

# Detect and model batch effects using SVA and limma
setwd(<>)

# Load GSE15573 data set using GEOquery, here I loaded saved 
# data object from my system
load("datasets/GSE15573/GSE15573_ESet.RData")
show(eset)

# Examine metadata
head(pData(eset))

# Select clinico-pathological variables from sample metadata
sample_data <- pData(eset)
sample_data <- data.frame( Trait = sample_data$characteristics_ch1,
                           Gender = sample_data$characteristics_ch1.1,
                           Age = sample_data$characteristics_ch1.2)

# Modify the variable definitions to make them more palatable
sample_data$Trait[sample_data$Trait == "status: Rheumatoid Arthritis Patient"] <- "RA"
sample_data$Trait[sample_data$Trait == "status: Control"] <- "Control"
sample_data$Gender <- gsub("gender: " , "", sample_data$Gender)
sample_data$Age <- gsub("age: " , "", sample_data$Age)
class(sample_data$Age) # age is still a character, we need to turn it into numbers
sample_data$Age <- as.numeric(sample_data$Age)
class(sample_data$Age)

# Convert age to categorical variable
age_group<-cut(sample_data$Age, 
               breaks=c(44, 58, 81), right = T)
age_group
sample_data$age_group <- age_group
table(sample_data$Trait, sample_data$age_group)

# Add modified variables to pheno data object
pData(eset)$Trait <- sample_data$Trait
pData(eset)$Gender <- sample_data$Gender
pData(eset)$Age <- sample_data$Age
pData(eset)$Age_group <- sample_data$age_group
head(pData(eset))

# Use SVA to identify latent sources of variance
pheno <- pData(eset)
head(pheno)

# Create null model matrix with adjustment variables only
mod0 <- model.matrix(~as.factor(Gender) + as.factor(Age_group), 
                     data = pheno)

# Create full model matrix that includes variable of interest and
# adjustment variables
mod <- model.matrix(~as.factor(Trait) + as.factor(Gender) + 
                      as.factor(Age_group), data = pheno)

# Detect latent batch effects
edata <- exprs(eset)
edata <- normalizeVSN(edata)
svobj <- sva(edata, mod, mod0)
head(svobj$sv)

# Use limma to correct for hidden batch effects detected with SVA
# Specify design matrix for the variables we want to preserve
mod <- model.matrix(~as.factor(Trait) + as.factor(Gender) + 
                      as.factor(Age_group), data = pheno)

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

# Cluster genes and samples, plot heatmap
pdf("GSE15573_heatmap_before_correction.pdf", width = 7, height = 8)
pheatmap(top500, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = sample_data[,-3], 
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
rownames(sample_data) <- colnames(edata_cor)

# Cluster genes and samples, plot heatmap
# Note that the samples are now clearly clustered according to the Trait
pdf("GSE15573_heatmap_after_correction.pdf", width = 7, height = 8)
pheatmap(top500, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = sample_data[,-3], 
         show_rownames = F)
dev.off()

# -------------------------- PCA ---------------------------------------- #
# Plot PCA after batch correction
p <- pca(edata_cor, metadata = sample_data[,-3], removeVar = 0.9)

# Create biplot and draw stat ellipses at 95% CI around groups
biplot(p,
       colby = 'Trait', colkey = c('Control' = 'forestgreen', 'RA' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)
ggsave("GSE15573_PCA_biplot_after_correct.pdf", device = "pdf", units = "in", 
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
ggsave("GSE15573_pairs_plot_after_correct.pdf", device = "pdf", units = "in", 
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
ggsave("GSE15573_loadings_plot_corrected.pdf", device = "pdf", units = "in", 
       width = 6, height = 7)

# Create eigencorrelation plot
pdf("GSE15573_eigencor_plot_correct.pdf", width = 7, height = 7) 
eigencorplot(p,
             metavars = c('Trait','Gender','age_group'))
dev.off()

# Identify optimal number of PCs
horn <- parallelPCA(expr) # Horn's method
horn$n # 7

elbow <- findElbowPoint(p$variance)
elbow # 6

# Create screeplot
screeplot(p, components = getComponents(p, 1:10),
          vline = c(horn$n, elbow)) + 
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))
ggsave("GSE15573_screeplot.pdf", device = "pdf", units = "in", 
       width = 7, height = 7)

# -------------------- Detect DEGs after correction ----------------------- #
# Create expression set object with VSN normalized data
# Note, that we should not use SVA corrected data in differential expression
# analysis. Instead, we should supply latent variables detected with SVA
# as covariates in model matrix
expr_vsn <- normalizeVSN(exprs(eset))
eset_vsn <- ExpressionSet(assayData = expr_vsn, phenoData = phenoData(eset), 
                           featureData = featureData(eset))
pheno <- pData(eset_vsn)
surrogate <- svobj$sv
colnames(surrogate) <- c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6")
pheno <- cbind(pData(eset_vsn), surrogate) # Add surrogate variables to pheno data 
head(pheno)
# Create model matrix with adjustment and surrogate variables
design <- model.matrix(~Gender + Age_group + Trait + SV1 +
                       SV2 + SV3 + SV4 + SV5 + SV6, 
                       data = pheno)
head(design)

# Fit the model
fit <- lmFit(eset_vsn, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "TraitRA", number = dim(fit$genes)[1], 
                adjust.method = "BH") 
head(res)
write.csv(res, file = "GSE15573_RA_vs_Controls_corrected.csv")

# How many genes are significant (adj.p.value < 0.05 and absolute 
# fold change > 1.5)
sum(res$adj.P.Val < 0.05) # 2720 probes
# Save the results for significant probes
write.csv(res[res$adj.P.Val < 0.05,], 
          file = "GSE15573_RA_vs_Controls_sig_only_corrected.csv")

# Plot a heatmap of differentially expressed probes
sig_probes <- res[which(res$adj.P.Val < 0.05),]$ID
expr <- exprs(eset_cor)
sig_exp <- expr[rownames(expr) %in% sig_probes,]

annot <- data.frame(Trait = pData(eset_cor)$Trait,
                    Gender = pData(eset_cor)$Gender,
                    Age_group = pData(eset_cor)$Age_group)
rownames(annot) <- rownames(pData(eset_cor))
head(annot)

# Clustering and heatmap for significant probes (adj. p-values < 0.05)
pdf("GSE15573_heatmap_significant_corrected.pdf", width = 7, height = 6)
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
ggsave("GSE15573_MAplot_corrected.pdf", device = "pdf", width = 7, height = 5)

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
ggsave("GSE15573_volcano_plot_corrected.pdf", device = "pdf", 
       width = 5, height = 7) 


