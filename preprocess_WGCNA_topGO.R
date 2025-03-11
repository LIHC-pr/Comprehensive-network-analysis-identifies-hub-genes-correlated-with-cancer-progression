library(UCSCXenaTools);
library(data.table);
library(R.utils);
library(dplyr);

setwd("E:/Dr Sharifi class project")

#Generate an object tracking all data sets from Xena Data Hubs.
data(XenaData);
write.csv(XenaData, "00_tblXenaHubInfo.csv")

#step5:Select then download target data sets from Xena Data Hubs
#[Step 5-a] Target = RSEM expected counts provided by the UCSC Toil Recompute Compendium
GeneExpectedCnt_toil = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count");
XenaQuery(GeneExpectedCnt_toil) %>%
  XenaDownload(destdir = "E:/Dr Sharifi class project")

#[Step 5-b] Target = TCGA Clinical data.
# Selecting the TCGA-LIHC cohort
paraCohort = "TCGA Liver Cancer";

# Selecting the appropriate dataset for TCGA-LIHC
# You can specify the dataset name related to BRCA. For example:
# Clinical matrix for BRCA could be something like "BRCA_clinicalMatrix" 
# (you'll need to confirm the exact dataset name from the Xena hub)
paraDatasets = "TCGA.LIHC.sampleMap/LIHC_clinicalMatrix";

# Generate the query for the TCGA-LIHC data
Clin_TCGA = XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterCohorts = paraCohort) %>%
  XenaFilter(filterDatasets = paraDatasets);

# Download the data to the current directory
XenaQuery(Clin_TCGA) %>%
  XenaDownload(destdir = "./")

#[Step 5-c] Target = TCGA Survival data.
Surv_TCGA = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%  
  XenaFilter(filterDatasets = "TCGA_survival_data");  
XenaQuery(Surv_TCGA) %>%
  XenaDownload(destdir = "./") 

#[Step 5-d] Target = GTEx Phenotype data.
Pheno_GTEx = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype");
XenaQuery(Pheno_GTEx) %>%
  XenaDownload(destdir = "./")


#[Step 6-a] Retrieve IDs for GTEx normal samples of desired tissue type(s).
filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz");
names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01));

paraStudy = "GTEX"; #Setting "GTEx" as the study of interest.
paraPrimarySiteGTEx = "Liver"; #Setting "Breast" as the primary site of interest.
paraPrimaryTissueGTEx = "Liver"; #Breast" as the primary tissue of interest.

filterGTEx02 = subset(filterGTEx01,
                      study == paraStudy &
                        primarysite == paraPrimarySiteGTEx &
                        grepl(paraPrimaryTissueGTEx, filterGTEx01$`primary disease or tissue`))

#[Step 6-b] Retrieve IDs for TCGA primary tumor samples of desired histological type(s).
filterTCGA01 = fread(paraDatasets);
names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01));

paraSampleType = "Primary Tumor"; #Setting "Primary Tumor" as the sample type of interest.
paraPrimarySiteTCGA = "Liver"; #Setting "Breast" as the primary site of interest.
paraHistologicalType = "Hepatocellular Carcinoma"; #Setting "Infiltrating Ductal Carcinoma" as the histological type of interest.

filterTCGA02 = subset(filterTCGA01,
                      sampletype == paraSampleType &
                        primarysite == paraPrimarySiteTCGA &
                        grepl(paraHistologicalType, filterTCGA01$histologicaltype))

#[Step 6-c] Merge GTEx and TCGA sample lists. Then pull expression profiles from .gz file by IDs on merged sample list.
filterExpr = c(filterGTEx02$sample, filterTCGA02$sampleID, "sample");

ExprSubsetBySamp = fread("TcgaTargetGtex_gene_expected_count.gz",
                         select = filterExpr)                         

#[Main Text: Step 7] Subset expression data to include only protein coding genes.
probemap = fread("zz_gencode.v23.annotation.csv", select = c(1, 2));
exprALL = merge(probemap, ExprSubsetBySamp, by.x = "id", by.y = "sample");
genesPC = fread("zz_gene.protein.coding.csv");
exprPC = subset(exprALL, gene %in% genesPC$Gene_Symbol);

#Remove duplicate gene symbols.
exprFinal = exprPC[!(duplicated(exprPC$gene) |
                       duplicated(exprPC$gene, fromLast = TRUE)), ]

#[Main Text: Step 8] Save the expression profile data frame for downstream analyses.
write.csv(exprFinal, "00_ExpectedCnt.csv")

#Print out a list of all available variables for the Liver Cancer cohort of TCGA.
names(filterTCGA02);

#Keep variable "pathologicstage".
varClinKeep = c("sampleID", "pathologicstage");
clinDF01 = as.data.frame(do.call(cbind, filterTCGA02));
clinFinal = clinDF01[varClinKeep];

#Identify observations/samples with no values assigned to the kept variables.
colSums(clinFinal == "");

colSums(is.na(clinFinal));

#Replace "no values" with "NA"
NA -> clinFinal[clinFinal == ""];
colSums(is.na(clinFinal));

#For the variable "pathologicstage", check count of YES/NO.
table(clinFinal$pathologicstage);

clinFinal <- clinFinal[clinFinal$pathologicstage != "[Discrepancy]", ]
table(clinFinal$pathologicstage);

#change the ordinal variable "pathologicstage" to numeric values.
clinFinal$pathologicstage[grepl("^Stage IV", clinFinal$pathologicstage)] <- 4;
clinFinal$pathologicstage[grepl("^Stage III", clinFinal$pathologicstage)] <- 3;
clinFinal$pathologicstage[grepl("^Stage II", clinFinal$pathologicstage)] <- 2;
clinFinal$pathologicstage[grepl("^Stage I", clinFinal$pathologicstage)] <- 1

#[Main Text: Step 10] Save the clinical data data frame for downstream analyses.
write.csv(clinFinal, "00_ClinTraits.csv")

table(clinFinal$pathologicstage);

#[Main Text: Step 11] Load R packages.
library(dplyr);
library(limma);
library(edgeR)

#[Main Text: Step 12] Back transformation of log-transformed expected count.
exprFinal = read.csv("00_ExpectedCnt.csv");
exprBT = exprFinal[, -c(1:3)];
rownames(exprBT) = exprFinal$gene;
exprBT = round(((2^exprBT)-1), 0);
write.csv(exprBT, "01_ExpectedCntBT.csv")

#[Main Text: Step 13] Convert count data to DGEList object.
expLIMMA = exprBT

x = DGEList(expLIMMA)

#[Main Text: Step 14] Group samples by condition (i.e., TCGA tumor or GTEx normal).
snames = colnames(x); 
group = substr(snames, 1, 4); #Sets up level information for samples.
x$samples$group = group #Assigns samples to appropriate group.

#[Main Text: Step 15] Compute counts per million.
cpm = cpm(x);
lcpm = cpm(x, log = TRUE);
L = mean(x$samples$lib.size) * 1e-6;
L; #Displays average library size in millions.
M = median(x$samples$lib.size) * 1e-6;
M; #Displays median library size in millions.
table(x$samples$group) #Returns number of samples per group.

#[Main Text: Step 16] Remove genes that are lowly expressed.
keep.exprs = filterByExpr(x, group = group);
x = x[keep.exprs, , keep.lib.sizes = FALSE];
dim(x) #Returns number of genes and samples retained.

#Generate density plot of log-CPM values for QC.
lcpm.cutoff = log2(10/M + 2/L);
nsamples = ncol(x);
par(mfrow=c(1, 2));
plot(density(lcpm[, 1]), lwd = 2, ylim = c(0, 0.26), las = 2, main = "", xlab = "");
title(main = "A. Expected Count", xlab = "Log-CPM");
abline(v = lcpm.cutoff, lty = 3);
for (i in 2:nsamples){
  den = density(lcpm[, i])
  lines(den$x, den$y, lwd = 2)
}
lcpm = cpm(x, log = TRUE);
plot(density(lcpm[, 1]), lwd = 2, ylim = c(0, 0.26), las = 2, main = "", xlab = "");
title(main = "B. Filtered Expected Count", xlab = "Log-CPM");
abline(v = lcpm.cutoff, lty = 3);
for (i in 2:nsamples){
  den = density(lcpm[, i])  
  lines(den$x, den$y, lwd = 2)  
}  

#[Main Text: Step 17] Compute scaling factors to convert observed library sizes into effective library sizes.
x = calcNormFactors(x, method = "upperquartile")

head(x$samples$norm.factors) #Prints out example normalization factors.

#[Main Text: Step 18] Generate design matrix.
design = model.matrix(~0 + group)

#[Main Text: Step 19] Set up contrast for comparison.
colnames(design) = gsub("group", "", colnames(design));
contr.matrix = makeContrasts(TCGAvsGTEX = TCGA - GTEX,
                             levels = colnames(design))

#[Main Text: Step 20] Transform RNA-Seq data for linear modeling.
v = voom(x, design, plot = TRUE)

#[Main Text: Step 21] Fit a linear model using weighted least squares for each gene.
vfit = lmFit(v, design)
vfit = contrasts.fit(vfit, contrasts = contr.matrix) #Returns the logFC between groups (TCGA vs. GTEx) via contra sts of the fitted linear models.

#[Main Text: Step 22] Perform empirical Bayes smoothing of standard errors.
efit = eBayes(vfit);
plotSA(efit, main="Final model: Mean-variance trend")

#[Main Text: Step 23] Examine the number of DEGs.
summary(decideTests(efit))

tfit = treat(vfit, lfc = 0.58);
summary(decideTests(tfit)) #Examines the number of DEGs with the added condition of minimum log-FC of 0.58.

#[Main Text: Step 24] Write out the DEG data frame for downstream analyses
DEGsTreat = topTreat(tfit, n = Inf);
write.csv(DEGsTreat, "01_DEGsTreat.csv")

#[Main Text: Step 25] Write out voom normalized expression profile for downstream WGCNA.
voomExpr = v$E;
write.csv(voomExpr, "01_voomExpr.csv")

#[Main Text: Step 26] Load R packages.
library(WGCNA);
library(dplyr)

#[Main Text: Step 27] For WGCNA, use only TCGA expression data.
options(stringsAsFactors = FALSE);
dataExpr0 = read.csv("01_voomExpr.csv");
head(dataExpr0[1:6]);

dataExpr1 = dataExpr0 %>% select(starts_with("TCGA"));
head(dataExpr1[1:6]);

rownames(dataExpr1) = dataExpr0$X;
head(dataExpr1[1:6])

#[Main Text: Step 28] Exclude genes with overall low expression.
#Add a new variable "Count" that counts the number of samples that have normalized expression value <0 for gene x.
dataExpr1$Count = rowSums(dataExpr1 < 0);
table(dataExpr1$Count);

#Keep genes that have "Count" = 0 (i.e., Keep genes that have non-negative normalized expression value across ALL samples).
dataExpr2 = dataExpr1 %>% filter(Count == 0);
#Remove the "Count" variable from the gene expression matrix.
dataExpr2$Count = NULL

#Transpose frame.
dataExpr3 = as.data.frame(t(dataExpr2))

#Double check for genes and samples with missing values. All genes and samples should be good (i.e., none needs to be removed). the output should be TRUE.
gsg = goodSamplesGenes(dataExpr3, verbose = 3);
gsg$allOK;

#[Main Text: Step 29] Cluster samples base on their Euclidean distance to see if there are any sample outliers.
sampleTree = hclust(dist(dataExpr3), method = "average");
plot(sampleTree,
     main = "Sample Clustering to Detect Sample Outliers",
     sub = "",
     xlab = "",
     cex = 0.6)

par(mfrow=c(1, 2));
byHist = hist(sampleTree$height,
              main = "Histogram of Height",
              xlab = "Height")

#height cutoff to remove outlier data is chosen based on dendrogram and histogram. best cutoff for this dataset is 90.
plot(sampleTree,
     main = "Sample Clustering to Detect Sample Outliers",
     sub = "",
     xlab = "",
     cex = 0.6);
abline(h = 100, col = "red"); #Cut tree at height 101.

clust = cutreeStatic(sampleTree,
                     cutHeight = 100,
                     minSize = 10);
table(clust) 

keepSamples = (clust == 1);
dataExpr4 = dataExpr3[keepSamples, ];
nGenes = ncol(dataExpr4);
nGenes; #Displays number of genes retained after removal of lowly-expressed genes.

nSamples = nrow(dataExpr4);
nSamples #Displays number of samples retained after removal of sample outliers.

#Re-cluster kept samples to inspect distribution
sampleTree2 = hclust(dist(dataExpr4), method = "average");
plot(sampleTree2,
     main = "Sample Clustering after Removal of Sample Outliers",
     sub = "",
     xlab = "",
     cex = 0.6);


dataExpr = dataExpr4

#[Main Text: Step 30] Prepare clinical trait data for WGCNA.
dataTraitALL = read.csv("00_ClinTraits.csv");
head(dataTraitALL);
dataTraitALL$sampleID = gsub("-", ".", dataTraitALL$sampleID);
head(dataTraitALL)

#Form a trait data frame analogous to the expression data frame.
exprRows = rownames(dataExpr);
traitRows = match(exprRows, dataTraitALL$sampleID);
dataTrait = dataTraitALL[traitRows, -1];
head(dataTrait);

rownames(dataTrait) = dataTrait$sampleID;
head(dataTrait);

dataTrait = select(dataTrait, -sampleID);
head(dataTrait)

#Visualize how the trait data relates to sample clustering.
sampleTree3 = hclust(dist(dataExpr), method = "average");
traitColors = numbers2colors(dataTrait, signed = FALSE)

#In the plot below, trait was converted to a color representation of absolute data value. Where white means low, red means high, and gray means missing entry.
plotDendroAndColors(sampleTree3,
                    traitColors,
                    groupLabels = names(dataTrait),
                    main = "Sample Dendrogram in Relation to Traits",
                    autoColorHeight = FALSE,
                    colorHeight = 0.1)

#[Save point] Backup TCGA gene expression and traits .Rdata.
save(dataExpr, dataTrait, file = "02_WGCNAdataInput.RData")

#[Main Text: Step 31] Analyze scale free topology for a set of soft-thresholding power to select appropriate soft thresholding power for network construction.
powers = c(c(1:10), seq(from = 12, to = 20, by = 2));
sft = pickSoftThreshold(dataExpr,
                        powerVector = powers,
                        networkType = "signed",
                        verbose = 5,
                        blockSize = 20000)

#[Main Text: Step 32] Plot results returned by the function pickSoftThreshold .
#Visualize scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale Independence"));
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels = powers,
     col = "red");
abline(h = 0.80,
       col = "red")

#Visualize the distribution of k(connectivity) for the selected soft-thresholding power (beta = 9).
k = softConnectivity(dataExpr,
                     type = "signed",
                     power = 9,
                     blockSize = 20000,
                     verbose = 2);

sum(is.na(k));

par(mfrow=c(1, 2));
byHist = hist(k)

#Visualize scale-free topology when the selected soft-thresholding power, beta = 9.
scaleFreePlot(k, main = "Scale Free Plot (Pearson), sft=9\n")

#[Main Text: Step 33] Calculate network adjacency from the gene expression matrix.
softPower = 9;
adjacency = adjacency(dataExpr,
                      power = softPower,
                      type = "signed")
#[Main Text: Step 34] Calculate the topological overlap matrix (TOM) from the adjacency matrix.
TOM = TOMsimilarity(adjacency,
                    TOMType = "signed",
                    verbose = 5)
dissTOM = 1 - TOM

#[Main Text: Step 35] Hierarchical clustering and module assignment.
geneTree = hclust(as.dist(dissTOM),
                  method = "average");
plot(geneTree,
     xlab = "",
     sub = "",
     main = "Gene Clustering",
     labels = FALSE,
     hang = 0.04);

minModuleSize = 75;
dynamicMods = cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            method = "hybrid",
                            deepSplit = 4,
                            pamStage = TRUE,
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize,
                            verbose = 4);

table(dynamicMods);

dynamicColors = labels2colors(dynamicMods);
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")

#[Main Text: Step 36] Calculate module eigengenes (MEs, i.e., first principal component).
MEList = moduleEigengenes(dataExpr,
                          colors = dynamicColors)

#Extract matrix containing the sample eigen value for each module.
MEs = MEList$eigengenes

#[Save point] Backup MEs and module assignments .Rdata.
save(MEs, dynamicMods, dynamicColors, geneTree,
     file = "02_WGCNAmodAssigned.RData")

#[Main Text: Step 37] Relate sample's eigen values to sample's clinical trait.
moduleTraitCor = cor(MEs, dataTrait, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor);
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dataTrait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = .65,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))

#[Main Text: Step 38] Relate gene expression levels to clinical trait(s). i.e., Gene Significance.
geneTraitSignificance = as.data.frame(cor(dataExpr, dataTrait$pathologicN, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = "GS_pathologicstage";
names(GSPvalue) = "p.GS_pathologicstage"

#[Main Text: Step 39] Relate gene expression levels to modules. i.e., Module Membership.
modNames = substring(names(MEs), 3);
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="")

#[Main Text: Step 40] Examine the gene significance (GS) vs. module membership (MM) relationship.
#Gene significance (GS) vs. module membership (MM) can be plotted for interesting modules. This is done to examine if genes that are highly associated with the trait of interest are also highly associated with their assigned module.
module = "red";
column = match(module, modNames);
moduleGenes = dynamicColors == module;
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "Module"),
                   ylab = "Gene Significance for pathologicstage",
                   main = paste ("Module Membership vs.Gene Significance\n"),
                   cex.main = 1,
                   cex.lab = 1,
                   cex.axis = 1,
                   col = module)

#[Main Text: Step 41] Annotate results from WGCNA.
annot = read.csv(file = "zz_gencode.v23.annotation.csv");
dim(annot);

head(annot);

probes = names(dataExpr);
probes2annot = match(probes, annot$gene);
sum(is.na(probes2annot)); #Checks for number of probes without annotation. Should be zero.

ExportPrep = data.frame(geneSymbol = probes,
                        geneSymbolCheck = annot$gene[probes2annot],
                        ENSG = annot$id[probes2annot],
                        moduleColor = dynamicColors,
                        geneTraitSignificance,
                        GSPvalue);
modOrder = order(-(cor(MEs, dataTrait$pathologicstage, use = "p")));
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(ExportPrep)
  ExportPrep = data.frame(ExportPrep, geneModuleMembership[, modOrder[mod]],
                          MMPvalue[, modOrder[mod]]);
  names(ExportPrep) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                        paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(ExportPrep$moduleColor, -abs(ExportPrep$GS_pathologicstage));
ExportFinal = ExportPrep[geneOrder, ]

#[Main Text: Step 42] Save annotated WGCNA results for downstream GO analysis.
write.csv(ExportFinal, file = "02_GSandMM.csv")

#Related to Gene set enrichment analysis with topGO Steps 43 - 51
#[Main Text: Step 43] Load R packages
library(data.table);
library(grex);
library(biomaRt);
library(topGO);
library(dplyr);
library(ggplot2)

#[Main Text: Step 44] Annotate results from WGCNA with Entrez IDs
annot = fread("02_GSandMM.csv",
              select = c("geneSymbol", "ENSG", "moduleColor"));
annot$ensembl = substr(annot$ENSG, 1, 15);
ensembl_id = annot$ensembl;
all_id = grex(ensembl_id);
annotComplete = merge(annot, all_id,
                      by.x = "ensembl",
                      by.y = "ensembl_id");
write.csv(annotComplete, "03_moduleGeneAnnotated.csv")

#[Main Text: Step 45] For GO annotation, connect to the ENSEMBL_MART_ENSEMBL BioMart databasae to query GO IDs for the input entrezgene_id .
genes_bg = annotComplete$entrez_id;
length(genes_bg);

tot_background = length(genes_bg);
db = useMart("ENSEMBL_MART_ENSEMBL",
             dataset = "hsapiens_gene_ensembl",
             host = "www.ensembl.org");
go_ids = getBM(attributes = c("go_id",
                              "entrezgene_id",
                              "namespace_1003"),
               filters = "entrezgene_id",
               values = genes_bg,
               mart = db)

#[Main Text: Step 46] Define (WGCNA) module of interest and set up named factors for genes of interest
modInt = as.factor(annotComplete$moduleColor);
annotSplit = split(annotComplete, modInt);
candidate_list = annotSplit$red$entrez_id;
length(candidate_list);

tot_candidate = length(candidate_list);
keep = candidate_list %in% go_ids[, 2];
keep = which(keep == TRUE);
candidate_list = candidate_list[keep];
geneList = factor(as.integer(genes_bg %in% candidate_list));
names(geneList) = genes_bg

#[Main Text: Step 47] Build topGO data object then run GO analysis
gene2GO = unstack(go_ids[, c(1, 2)]);
GOdata = new("topGOdata",
             ontology = c("BP"),
             allGenes = geneList,
             annot = annFUN.gene2GO,
             gene2GO = gene2GO);
allGO = usedGO(GOdata)

#[Main Text: Step 48] Test significance of enriched GO terms.
sigTest = runTest(GOdata,
                  algorithm = "elim",
                  statistic = "fisher")

#[Main Text: Step 49] Generate a summary table of topGO results.
all_res = GenTable(GOdata,
                   weightFisher = sigTest,
                   orderBy = "weightFisher",
                   topNodes = length(allGO))

#[Main Text: Step 50] Calculate odds ratios.
all_res$OR = log2((all_res$Significant/tot_candidate)/(all_res$Annotated/tot_background))
write.csv(all_res, "03_moduleGOAnnotated.csv")

#[Main Text: Step 51] Generate a summary figure of topGO results.
#For this example, threshold was set at weightFisher < 0.05 and annotated.backgound was set at Annotated >= 30 . As well, only the top 10 most significant GO terms were selected for presentation.
#This sets the threshold p < 0.05 and #.annotated.background.genes >= 30.
GO_bar = all_res %>%
  filter(weightFisher < 0.05) %>%
  filter(Annotated >= 30);
GO_bar = GO_bar %>%
  dplyr::select(Term, weightFisher, OR);
sapply(GO_bar, class);

GO_bar = transform(GO_bar, weightFisher = as.numeric(weightFisher));
sapply(GO_bar, class);

#This selects the top 10 most significant (i.e., order by ascending p-values) GO terms for presentation.
GO_bar = GO_bar %>%
  arrange(weightFisher)

CapStr = function(y) {
  c = strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1, 1)), substring(c, 2),
        sep = "", collapse = " ")
}
GO_bar$Term = sapply(GO_bar$Term, CapStr);
#This selects the top 10 most significant (i.e., order by ascending p-values) GO terms for presentation.
GO_bar = GO_bar %>%
  arrange(weightFisher) %>%
  head(20);
CapStr = function(y) {
  c = strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1, 1)), substring(c, 2),
        sep = "", collapse = " ")
}
GO_bar$Term = sapply(GO_bar$Term, CapStr);
p = ggplot(GO_bar,
           aes(x = reorder(Term, OR), y = OR, fill = weightFisher)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  coord_flip() +
  scale_fill_gradient(low="#feff2b",high="#fe0100") +
  ylim(0, 3) +
  labs(title = ~underline("Enriched GO Biological Processes"),
       x = NULL,
       y = "Odds Ratio",
       fill = "p.value") +
  theme_bw() +
  theme(plot.title.position = "plot") +
  theme(plot.title = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold")) +
  theme(legend.position = "right");
p
