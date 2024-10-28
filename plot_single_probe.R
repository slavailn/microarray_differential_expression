# Plot the expression of single probe as boxplot overlaid with a stripchart
# The data is represented as Expression Set object
# The treatment groups are "control" and "dose5Gy", replace these values depending on the dataset
library(limma)
idx_ct <- which(pData(normData)$exposure == "control")
idx_tr <- which(pData(normData)$exposure == "dose5Gy") 
idx <- c(idx_ct, idx_tr)
exp <- exprs(normData)[,idx]
exp <- exp[which(rownames(exp) == "1421037_at"),]
exp <- data.frame(row.names = names(exp), 
                  expression = exp,
                  exposure = pData(normData)$exposure[idx])
boxplot(expression ~ exposure, data = exp, lwd = 2, ylab = 'expression',
        main = "1421037_at")
stripchart(expression ~ exposure, vertical = TRUE, data = exp, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
