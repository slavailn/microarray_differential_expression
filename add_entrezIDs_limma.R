library("hgu133plus2.db")

setwd(<path/to/wd>)

# Annotate limma results with EntrezIDs for GPL570 platform
# with hgu133plus2.db package
x <- hgu133plus2ENTREZID
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
xx <- as.data.frame(unlist(xx))
names(xx) <- "EntrezID"

addAnnot <- function(resFile, annotDF, outName) 
{
    # Get results
    res <- read.csv(resFile, header = T, 
                row.names = 1)
    res_annot <- merge(res, xx, by = 0) 
    # Sort by p-value
    res_annot <- res_annot[order(res_annot$P.Value),]
    # Remove duplicated EntrezIDs, keep only the probe with the lowest p-value
    res_annot <- res_annot[!duplicated(res_annot$EntrezID),]
    write.csv(res_annot, outName, row.names = F)
}

addAnnot(resFile = "results/res.csv",
         annotDF = xx, outName = "res.annot.csv")




