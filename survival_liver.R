getwd()

#[Main Text: Step 52] Load R packages.
library(dplyr);
library(data.table);
library(RegParallel);
library(survminer)

#[Main Text: Step 53] Subset the voom normalized expression data to include only TCGA sample IDs and only genes within the module of interest.
dataExpr = read.csv("01_voomExpr.csv");
rownames(dataExpr) = dataExpr$X;
dataExpr = dataExpr %>% dplyr::select(starts_with("TCGA"));
head(dataExpr[1:4]);

genes_toKeep = read.csv("02_GSandMM.csv");
head(genes_toKeep[1:5]);

genes_toKeep_list = genes_toKeep %>% filter(moduleColor == "red");
genes_toKeep_list = genes_toKeep_list$geneSymbol;
dataExpr = dataExpr[genes_toKeep_list, ];
head(dataExpr[1:4])

#[Main Text: Step 54] Transform expression data to Z-score.
#Also, set Z-score cut-offs for high and low expression. Z-score of (+/-)0.674 translates to 25/75 percentile.
dataExpr = as.data.frame(t(dataExpr));
dataExpr = as.data.frame(scale(dataExpr));
head(dataExpr[1:4]);

highExpr = 0.674;
lowExpr = -0.674;
dataExpr = as.data.frame(ifelse(dataExpr <= lowExpr, 1,
                                ifelse(dataExpr >= highExpr, 2, NA)));
dataExpr[] = lapply(dataExpr, factor);
dataExpr$SampleID = rownames(dataExpr)

#[Main Text: Step 55] Create an object that contains both survival and Z-score data.
dataSurv = as.data.frame(fread("TCGA_survival_data"));
dataSurv$sample = gsub("-", ".", dataSurv$sample);
head(dataSurv);

dataJoined = merge(dataSurv, dataExpr, by.x = "sample", by.y = "SampleID");

#Generate a histogram that summarize the range and distribution of available survival data.
temp4Hist = (dataJoined$OS.time)/365;
hist(temp4Hist,
     breaks = seq(from = 0, to = 26, by = 1),
     xaxt='n',
     xlab = "Years",
     main = "Histogram of OS.time");
axis(side = 1,
     at = seq(0, 26, 1),
     labels = seq(0, 26, 1));

#Histogram tells us that data is available for examining three-, five-, and/or 10-year overall survival probability. Subset OS.time for examination of 10-year overall survival probability
dataJoined = subset(dataJoined, OS.time <= 3650)


#[Main Text: Step 56] Run RegParallel to test each gene independently via Cox regression.
res = RegParallel(data = dataJoined,
                  formula = 'Surv(OS.time, OS) ~ [*]',
                  FUN = function(formula, data) coxph(formula = formula,
                                                      data = data,
                                                      ties = 'breslow',
                                                      singular.ok = TRUE),
                  FUNtype = 'coxph',
                  variables = colnames(dataJoined)[10:ncol(dataJoined)],
                  blocksize = 151,
                  cores = 2,
                  nestedParallel = FALSE,
                  conflevel = 95);
res = res[!is.na(res$P), ];
res = res[order(res$LogRank, decreasing = FALSE)];
head(res);


write.csv(res, "04_Survival.csv")

#[Main Text: Step 57] Load in the summary of limma-voom differential expression analysis.
dataDEG = read.csv("01_DEGsTreat.csv");
head(dataDEG)

#[Main Text: Step 58] Create nodes and edges files for use in Cytoscape.
dataCytoscape = merge(res, dataDEG, by.x = "Variable", by.y = "X");
dataCytoscape = dataCytoscape %>% dplyr::select(Variable, LogRank, logFC, adj.P.Val);
head(dataCytoscape);

write.csv(dataCytoscape, "04_CytoscapeInput.csv")

#[Main Text: Step 59] Generate a short-list of genes that satisfy the conditions:
#a) Is associated with clinical trait of interest (e.g., lymphatic invasion), b) is differentially expressed, and c) has significant effect on survival.
shortList = subset(dataCytoscape,
                   LogRank < 0.05 & adj.P.Val < 0.05);
nrow(shortList) #Quick inspection on number of genes that made the short-list.

#[Main Text: Step 60] Generate Kaplan-Meier (KM) plots for shortlisted genes.
geneSL = shortList$Variable;
for(i in 1:length(geneSL)) {
  fit = survfit(as.formula(paste0("Surv(OS.time, OS) ~", geneSL[i])),
                data = dataJoined)
  print(
    ggsurvplot(fit,
               pval = TRUE,
               risk.table = TRUE,
               break.time.by = 365,
               ggtheme = theme_bw(),
               palette = c("blue", "red"),
               xlim = c(0, 3650),
               risk.table.y.text.col = TRUE,
               risk.table.y.text = FALSE,
               tables.theme = theme_cleantable(),
               tables.height = 0.15,
               xlab = "Days",
               font.x = c(12, "bold"),
               ylab = "Survival Probability",
               font.y = c(12, "bold"),
               legend.labs = c("Low-expression", "High-expression"),
               font.legend = c(12, "bold"),
               legend.title = geneSL[i])
  )
}
