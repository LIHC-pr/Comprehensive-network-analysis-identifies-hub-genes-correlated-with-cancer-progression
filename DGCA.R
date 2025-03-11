getwd()
setwd("E:/Dr Sharifi class project")

library(dplyr);
dataExpr = read.csv("01_voomExpr.csv");
head(dataExpr[1:6]);

rownames(dataExpr) = dataExpr$X;
dataExpr = dataExpr[, -1];
head(dataExpr[1:6]);

dataSample = names(dataExpr);
group = substr(dataSample, 1, 4);
dataGTEx = as.numeric(ifelse(group == "GTEX", 1, 0));
dataTCGA = as.numeric(ifelse(group == "TCGA", 1, 0));
dataMatrix = data.frame(dataGTEx, dataTCGA);
dataMatrix = as.matrix(dataMatrix);

library(DGCA);
nrow(dataExpr);

dataExpr = filterGenes(dataExpr,
                       filterTypes = "central",
                       keepRows = NULL,
                       filterCentralType = "median",
                       filterCentralPercentile = 0.25,
                       allGroups = FALSE,
                       design = NULL);
nrow(dataExpr);

dataExpr = filterGenes(dataExpr,
                       filterTypes = "dispersion",
                       keepRows = NULL,
                       filterDispersionType = "dispersion_index",
                       filterDispersionPercentile = 0.25,
                       allGroups = FALSE,
                       design = NULL);
nrow(dataExpr);


#Compute the differential correlation (Pearson) of ESRP2 vs. All.
ddcor_BDH2_PS = ddcorAll(inputMat = dataExpr,
                             design = dataMatrix,
                             inputMatB = NULL,
                             compare = c("dataGTEx", "dataTCGA"),
                             splitSet = "BDH2",
                             impute = FALSE,
                             corrType = "pearson",
                             nPairs = "all",
                             sortBy = "zScoreDiff",
                             adjust = "perm",
                             nPerms = 1000,
                             classify = TRUE,
                             sigThresh = 0.05,
                             corSigThresh = 0.05,
                             heatmapPlot = FALSE,
                             color_palette = NULL,
                             verbose = FALSE,
                             plotFdr = FALSE,
                             corr_cutoff = 0.99,
                             signType = "none",
                             getDCorAvg = FALSE,
                             oneSidedPVal = FALSE);

#GO enrichment analysis with DGCA/GOstats.
library(org.Hs.eg.db);
library(GOstats);
library(HGNChelper);
library(plotrix);

ddcorGO_BDH2_PS = ddcorGO(ddcor_BDH2_PS,
                              universe = rownames(dataExpr),
                              pval_gene_thresh = 0.05,
                              classes = FALSE,
                              geneNameCol = "Gene1",
                              pval_GO_cutoff = 1,
                              HGNC_clean = TRUE,
                              HGNC_switch = TRUE,
                              gene_ontology = "all",
                              adjusted = TRUE,
                              annotation = "org.Hs.eg.db",
                              conditional = TRUE,
                              calculateVariance = TRUE,
                              unique_genes = FALSE,
                              regcor = FALSE,
                              ddcor_find_significant = TRUE,
                              ddcorGO_res = NULL);

plotGOTwoGroups(dfList1 = ddcorGO_BDH2_PS$enrichment_significant_gain_of_correlation_genes,
                dfList2 = ddcorGO_BDH2_PS$enrichment_significant_loss_of_correlation_genes,
                nTerms = 5,
                minSize = 30,
                maxSize = 1000,
                labelsCol = "Ontology",
                adjustPVals = TRUE,
                plotrix_gap = 9,
                GOTermTypes = c("BP", "CC", "MF"),
                pValCutoff = 0.05,
                filterSignificant = TRUE,
                filterSigThresh = 0.01,
                labels = c("Gain in Correlation",
                           "GO Term",
                           "Loss in Correlation"),
                fill_zero_cats = FALSE)
