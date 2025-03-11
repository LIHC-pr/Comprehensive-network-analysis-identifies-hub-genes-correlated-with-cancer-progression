setwd("D:/Dr sharifi zarchi bioinf class project")
library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(umap)

series <- "GSE54236"
platform <- "GPL6480"
gset <- getGEO("GSE54236", GSEMatrix =TRUE,AnnotGPL=TRUE)
gr <- c(rep("tumor", 81), rep("non",80))
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)

# log2 transformation
ex[is.na(ex)] <- min(ex, na.rm = T) #data imputation 
#exprs(gset) <- log2(ex - min(ex)+1)

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000001111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111")
sml <- strsplit(gsms, split="")[[1]]

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("tumor","non"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- subset(tT, select=c("adj.P.Val","logFC","Gene.symbol"))
tT2<- subset(tT, abs(logFC)>1 & adj.P.Val<0.05)
up.genes<-subset(tT2, logFC>0)
down.genes<-subset(tT2, logFC<0)
write.table(up.genes$Gene.symbol, file=stdout(), row.names=F,col.names =F ,quote=F, sep="\t")
write.table(down.genes$Gene.symbol, file=stdout(), row.names=F,col.names =F ,quote=F, sep="\t")

#box plot:
ex <- exprs(gset)
exm <- melt(ex)
colnames(exm) <- c("Gene", "Sample", "Expression")
ggplot(exm, aes(Sample, Expression, fill = Sample)) +
  geom_boxplot(width = 0.9) +
  theme(axis.text.x = element_text(angle = 90,size = 12,  hjust = 0.5,vjust = 1,lineheight = 1.2),plot.margin = margin(5, 5, 5, 5), legend.text = element_text(size = 8),  # Smaller legend text
        legend.title = element_text(size = 10)) +
  scale_x_discrete(expand = c(0, 0)) +  # Remove axis padding
  labs(x = "Sample", y = "Expression", fill = "Sample") +
  guides(fill = guide_legend(ncol = 3))  # Legend in 3 columns
ggsave("box123.pdf", width = 50, height = 15, device = "pdf", dpi = 300, limitsize = FALSE)

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")


# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE54236", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")
#PCA
cnt <- normalizeBetweenArrays(exprs(gset), method = "quantile")
ex  <- log2(1+cnt)
pc <- prcomp(ex - rowMeans(ex))
pcr <- data.frame(pc$r, group = gs)

#ggplot(pcr , aes(pc1, pc2))

#pheatmap
pheatmap(cor(ex - rowMeans(ex)), labels_row = gs, labels_col = gs, color = bluered(256))

# Select differentially expressed genes
selected_genes <- rownames(tT2) 
heatmap_data <- ex[selected_genes, ]

# Plot Heatmap
pheatmap(
  heatmap_data,
  scale = "row",              
  cluster_rows = TRUE,     
  cluster_cols = TRUE,        
  show_rownames = TRUE,      
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "black", "red"))(5000)
)