library(limma)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
# This script is designed to detect differentially 
# expressed genes in Affymetrix array with limma assuming that the data
# is already normalized

# Detect and model batch effects using SVA and limma
setwd("<working_directory>")
load("study_data/GSE15645_affy_norm_outlier_removed.RData")

# --------------------------------------------------------------------# 
# Examine metadata

# Create model matrix to compare RA samples to controls, 
# while accounting for biological factors of no interest, i.e.
# age and gender
# The variables need to be converted to factors
pData(Affy_norm)$Trait <- factor(pData(Affy_norm)$Trait, 
                                 levels = c("control", "polyarticular JIA"))
# Note, that we set Control as a reference level
pData(Affy_norm)$Trait

# Convert other variables to factors
pData(Affy_norm)$Disease_state <- factor(pData(Affy_norm)$Disease_state)
pData(Affy_norm)$Disease_state

# Create design matrix
design <- model.matrix(~Trait, data=pData(Affy_norm))
head(design)

# Fit the model
fit <- lmFit(Affy_norm, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "Traitpolyarticular JIA", number = dim(fit$p.value)[1], 
                adjust.method = "BH") 
head(res)
write.csv(res, file = "GSE15645_RA_vs_Controls.csv")

# How many genes are significant (adj.p.value < 0.05 and absolute 
# fold change > 1.5)
sum(res$adj.P.Val < 0.05) # 8108 probes
# Save the results for significant probes
write.csv(res[res$adj.P.Val < 0.05,], 
          file = "GSE15645_RA_vs_Controls_sig_only.csv")


# Plot a heatmap of differentially expressed probes
sig_probes <- rownames(res[which(res$adj.P.Val < 0.05),])
expr <- exprs(Affy_norm)
sig_exp <- expr[rownames(expr) %in% sig_probes,]
colnames(sig_exp) <- gsub(".CEL", "", colnames(sig_exp))

annot <- data.frame(Trait = pData(Affy_norm)$Trait,
                    Disease_state = pData(Affy_norm)$Disease_state)
rownames(annot) <- rownames(pData(Affy_norm))
head(annot)

# Clustering and heatmap for significant probes (adj. p-values < 0.05)
pdf("GSE15645_heatmap_significant.pdf", width = 7, height = 6)
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
ggsave("GSE15645_MAplot.pdf", device = "pdf", width = 7, height = 5)

# Create volcano plot
res <- res[with(res, order(logFC)),]
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
ggsave("GSE15645_volcano_plot.pdf", device = "pdf", width = 5, height = 7)
