#load dataset
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
gse <- getGEO("GSE5281",GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
expression_data <- exprs(gse[[1]])
#preprocessing
cat("Missing values in expression data:", sum(is.na(expression_data)), "\n")
#before normalisation
png(file = file.path("boxplot_before_normalization.png"))
boxplot(expression_data, main = "Before Normalization", las = 2, outline = FALSE)
dev.off()
#after normalisation
expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")
png(file = file.path("boxplot_after_normalization.png"))
boxplot(expression_data, main = "After Normalization", las = 2, outline = FALSE)
dev.off()
# Filter low-expressed genes
gene_means <- rowMeans(expression_data)
threshold <- quantile(gene_means, probs = 0.25)
expression_data <- expression_data[gene_means > threshold, ]
cat("Dimensions after filtering:", dim(expression_data), "\n")

metadata <- pData(gse[[1]])
colnames(metadata)
head(metadata)
# Extract clean columns
meta <- data.frame(
  sample = rownames(metadata),
  disease = metadata$`Disease State:ch1`,
  region  = metadata$`Organ Region:ch1`,
  sex     = metadata$`Sex:ch1`,
  age     = metadata$`Age:ch1`
)

head(meta)
#data framing acording to Entorhinal Cortex (EC) with Control vs Alzheimer
meta <- data.frame(
  sample = rownames(metadata),
  disease = metadata$`Disease State:ch1`,
  region  = metadata$`Organ Region:ch1`
)
#filter out na missing file
sum(is.na(meta$disease))
sum(is.na(meta$region))
#viewing missing samples
meta[is.na(meta$region) | is.na(meta$disease), ]
#remove na 
meta_clean <- meta[!is.na(meta$disease) & !is.na(meta$region), ]
head(meta_clean)

idx <- grepl("Entorhinal Cortex", meta$region)
meta_ec <- meta[idx, ]
table(meta_ec$disease)

head(meta_ec)
expr_ec <- expression_data[, idx]
dim(expr_ec)
#idexing
idx <- grepl("Entorhinal Cortex", meta$region)
sum(idx)   # should be 21 (you already confirmed)
meta_ec <- meta[idx, ]
expr_ec <- expression_data[, idx]
#rowname = colname
rownames(meta_ec) <- meta_ec$sample
head(colnames(expr_ec))
head(meta_ec$sample)
table(meta_ec$disease)

# Rename safely
group <- factor(group,
                levels = c("normal", "Alzheimer's Disease"),
                labels = c("Control", "AD"))
levels(group)
#limma design
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

design
contrast <- makeContrasts(AD - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf)
results <- topTable(fit2, number = Inf)

head(results)
deg <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
dim(deg)
#differential gene analysis
library(ggplot2)
#valcano plot
results$significant <- results$adj.P.Val < 0.05 & abs(results$logFC) > 1

ggplot(results, aes(logFC, -log10(adj.P.Val), color = significant)) +
  geom_point() +
  theme_minimal()
#heatmap
library(pheatmap)

top50 <- rownames(deg)[1:50]
pheatmap(expr_ec[top50, ], scale = "row")
#get into strict analysis to choose biomarker
deg_strict <- results[
  results$adj.P.Val < 0.01 & abs(results$logFC) > 2,
]

dim(deg_strict)
#from this pick out top genes 
top_genes <- head(deg_strict[order(deg_strict$adj.P.Val), ], 10)
top_genes
#from top gene choose one gen
gene <- "205914_s_at"
boxplot(expr_ec[gene, ] ~ meta_ec$disease,
        main = gene,
        col = c("blue", "red"))
#convert probe id to gene symbol
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)

library(AnnotationDbi)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = rownames(deg_strict),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg_strict$symbol <- gene_symbols
head(deg_strict)
#roc analysis
install.packages("pROC")
library(pROC)
roc1 <- roc(meta_ec$disease, expr_ec["205914_s_at", ])

roc2 <- roc(meta_ec$disease, expr_ec["214333_x_at", ])
roc3 <- roc(meta_ec$disease, expr_ec["209364_at", ])

plot(roc1, col="red")
plot(roc2, col="blue", add=TRUE)
plot(roc3, col="green", add=TRUE)

legend("bottomright",
       legend=c("205914", "214333", "209364"),
       col=c("red","blue","green"),
       lwd=2)
#auc numericaly
auc(roc1)  # for 205914
auc(roc2)  # for 214333
auc(roc3)  # for 209364
coords(roc3, "best", ret=c("threshold","sensitivity","specificity"))

deg_strict[deg_strict$symbol == "BAD", ]
head(deg_strict$symbol)
