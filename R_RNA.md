---
title: "R_RNA"
output: html_document
---


by adding `fig.path = 'figures/'` we put all of the figures created when we knit this document into a directory called `figures`


# Differential Expression Testing

Read the docs: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

Installs:

```r
#bioClite(GSEAbase)
#bioClite("clusterProfiler")
#install.packages("devtools")
#install.packages("RColorBrewer")
#install.packages("pheatmap")
#devtools::install_github("karthik/wesanderson")
#biocLite("org.Sc.sgd.db")
#biocLite("GOstats")
#biocLite(edgeR)
#install.packages("treemap")
```


Load Libraries: 

```r
library(tximport)
library(DESeq2)
library(tidyverse)
library("GSEABase")
library(clusterProfiler)
library(RColorBrewer)
library(pheatmap)
library(wesanderson)
library(org.Sc.sgd.db)
library(GOstats)
library(edgeR)
library(treemap)
```

Import sample metadata: 

```r
# read in the file from url
samples <- read_csv("https://osf.io/cxp2w/download")
# look at the first 6 lines
samples
```

```
## # A tibble: 6 x 3
##   sample    quant_file                                condition
##   <chr>     <chr>                                     <chr>    
## 1 ERR458493 ~/quant/ERR458493.qc.fq.gz_quant/quant.sf wt       
## 2 ERR458494 ~/quant/ERR458494.qc.fq.gz_quant/quant.sf wt       
## 3 ERR458495 ~/quant/ERR458495.qc.fq.gz_quant/quant.sf wt       
## 4 ERR458500 ~/quant/ERR458500.qc.fq.gz_quant/quant.sf snf2     
## 5 ERR458501 ~/quant/ERR458501.qc.fq.gz_quant/quant.sf snf2     
## 6 ERR458502 ~/quant/ERR458502.qc.fq.gz_quant/quant.sf snf2
```

Import tx 2 gene file: 

```r
tx2gene_map <- read_tsv("https://osf.io/a75zm/download")
```

```
## Parsed with column specification:
## cols(
##   TXNAME = col_character(),
##   GENEID = col_character()
## )
```

```r
txi <- tximport(files = samples$quant_file, type = "salmon", tx2gene = tx2gene_map)
```

```
## reading in files with read_tsv
```

```
## 1 2 3 4 5 6 
## summarizing abundance
## summarizing counts
## summarizing length
```

```r
colnames(txi$counts) <- samples$sample
```

Make DESeq2 object: 

```r
dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = samples, 
                                design = ~condition)
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
## design formula are characters, converting to factors
```

```
## using counts and average transcript lengths from tximport
```

```r
dds$condition <- relevel(dds$condition, ref = "wt") # make wild-type the reference to which expression in treatment samples is compared to 
```

Run DESeq2: 

```r
dds <- DESeq(dds)
```

```
## estimating size factors
```

```
## using 'avgTxLength' from assays(dds), correcting for library size
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

Check out results: 

```r
res <- results(dds)
head(res)
```

```
## log2 fold change (MLE): condition snf2 vs wt 
## Wald test p-value: condition snf2 vs wt 
## DataFrame with 6 rows and 6 columns
##          baseMean log2FoldChange     lfcSE       stat       pvalue
##         <numeric>      <numeric> <numeric>  <numeric>    <numeric>
## ETS1-1 150.411381      0.1994166 0.1488952  1.3393084 1.804703e-01
## ETS2-1   0.000000             NA        NA         NA           NA
## HRA1     1.088842     -3.0083691 1.9214252 -1.5656967 1.174196e-01
## ICR1    28.871331      0.1815095 0.3479181  0.5217017 6.018780e-01
## IRT1    29.605385     -0.4930828 0.3335646 -1.4782227 1.393482e-01
## ITS1-1  23.867618      7.4605256 1.2123998  6.1535196 7.578202e-10
##                padj
##           <numeric>
## ETS1-1 3.563935e-01
## ETS2-1           NA
## HRA1             NA
## ICR1   7.596686e-01
## IRT1   2.972488e-01
## ITS1-1 9.351823e-09
```

Summarize results

```r
summary(res, alpha = 0.05) # default significance cut-off is 0.1, changing alpha to 0.05 changes the significance cut-off 
```

```
## 
## out of 6041 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)     : 598, 9.9% 
## LFC < 0 (down)   : 1049, 17% 
## outliers [1]     : 6, 0.099% 
## low counts [2]   : 235, 3.9% 
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

# Visualizing RNA-seq results 

## Normalization

**Count Data Transformations:** 
for ranking and visualizations (e.g. PCA plots and heatmaps)

**rlog**: "transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. The rlog transformation produces a similar variance stabilizing effect as varianceStabilizingTransformation, though rlog is more robust in the case when the size factors vary widely. The transformation is useful when checking for outliers or as input for machine learning techniques such as clustering or linear discriminant analysis." -- from function documentation 

This is computationally very time intensive. 


```r
rld <- rlog(dds, blind=TRUE)
head(assay(rld), 3)
```

```
##          ERR458493   ERR458494  ERR458495  ERR458500  ERR458501 ERR458502
## ETS1-1 7.194109821  7.20820460  7.1308260  7.3648267  7.2916848  7.189450
## ETS2-1 0.000000000  0.00000000  0.0000000  0.0000000  0.0000000  0.000000
## HRA1   0.004389821 -0.08089499 -0.1259611 -0.1080032 -0.1366897 -0.136613
```

** Variance stabilizing transformation (so much faster than rlog):**
"This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values). The transformation also normalizes with respect to library size. The rlog is less sensitive to size factors, which can be an issue when size factors vary widely. These transformations are useful when checking for outliers or as input for machine learning techniques such as clustering or linear discriminant analysis."" – from function documentation


```r
vsd <- vst(dds, blind = TRUE)
head(assay(vsd), 3)
```

```
##        ERR458493 ERR458494 ERR458495 ERR458500 ERR458501 ERR458502
## ETS1-1  7.414166  7.437301  7.308989  7.681222  7.568717  7.409575
## ETS2-1  3.842069  3.842069  3.842069  3.842069  3.842069  3.842069
## HRA1    4.629478  4.303053  3.842069  4.154975  3.842069  3.842069
```

## Ordination

rlog PCA: 

```r
data1 <- plotPCA(rld, returnData=TRUE)
data1$group<-gsub(" : ","_",as.character(data1$group))
percentVar1 <- round(100 * attr(data1, "percentVar"))

PCA<-ggplot(data1, aes(PC1, PC2, color = condition))+ theme_bw()+
  geom_point(size=9, alpha = 0.8) + scale_colour_manual(values = c("#44aabb","#bbbbbb"))+
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) +
  theme(text = element_text(size=20)) + ggtitle("rlog PCA")
PCA
```

![plot of chunk pca_rld](figures/pca_rld-1.png)

```r
#ggsave("figures/vsd_PCA.png", device="png") # to save the plot
```

variance stabilized PCA:

```r
data1 <- plotPCA(vsd, returnData=TRUE)
data1$group<-gsub(" : ","_",as.character(data1$group))
percentVar1 <- round(100 * attr(data1, "percentVar"))

PCA<-ggplot(data1, aes(PC1, PC2, color = condition))+ theme_bw()+
  geom_point(size=9, alpha = 0.8) + scale_colour_manual(values = c("#44aabb","#bbbbbb"))+
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) +
  theme(text = element_text(size=20)) + ggtitle("vst PCA")
PCA
```

![plot of chunk pca_vst](figures/pca_vst-1.png)

```r
#ggsave("figures/vsd_PCA.png", device="png") # to save the plot
```

## HeatMaps

rlog HeatMap:

```r
df <- as.data.frame(colData(rld)[,c("condition", "sample")])

mat_colors1<-list(sample = brewer.pal(12, "Paired")[0:6])
names(mat_colors1$sample)<- df$sample

mat_colors <- list(condition = brewer.pal(12, "Paired")[7:8])
names(mat_colors$condition) <- c("wt", "snf2")

genes <- order(res$padj)[1:1000]

 pheatmap(assay(rld)[genes, ], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, annotation_colors = c(mat_colors1, mat_colors), fontsize = 12)
```

![plot of chunk heatmap_rld](figures/heatmap_rld-1.png)

variance stabilized HeatMap: 

```r
df <- as.data.frame(colData(vsd)[,c("condition", "sample")])

pheatmap(assay(vsd)[genes, ], cluster_rows=TRUE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df, annotation_colors = c(mat_colors1, mat_colors), fontsize = 12)
```

![plot of chunk heatmap_vst](figures/heatmap_vst-1.png)

Another option for heat maps: 
plot the difference from the mean normalized count across samples 
(and optionally change default colors)

With Rlog transformed data:

```r
library(wesanderson)
pal <- wes_palette(name = "Zissou1", n=2000 , type= "continuous")

mat_colors1<-list(sample = wes_palette("IsleofDogs1", 6))
names(mat_colors1$sample)<- df$sample

mat_colors <- list(condition = wes_palette("Cavalcanti1")[4:5])
names(mat_colors$condition) <- c("wt", "snf2")

mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(rld)[,c("condition", "sample")])

pheatmap(mat,  cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, annotation_colors = c(mat_colors1, mat_colors), fontsize = 12, color = pal)
```

![plot of chunk heatmap_rld_meandiff](figures/heatmap_rld_meandiff-1.png)

Same but with variance stabilizing function:

```r
mat <- assay(vsd)[genes, ]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(vsd)[,c("condition", "sample")])

pheatmap(mat,  cluster_rows=TRUE, show_rownames=FALSE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df, annotation_colors = c(mat_colors1, mat_colors), fontsize = 12, color = pal)
```

![plot of chunk heatmap_vst_meandiff](figures/heatmap_vst_meandiff-1.png)


Heatmap of sample-to-sample distances

```r
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

![plot of chunk heatmap_sampledistance](figures/heatmap_sampledistance-1.png)

# Gene Set Enrichment Testing 
If you remember, we had  598 significantly upregulated genes and 1049 significantly down regulated genes in this data set (this is pretty typical). That is a lot to try to make sense of. If you know you are interested in a specific gene or a specific pathway, you can look for that in your data, but if you are trying to figure out what is generally different betwene treatments, it helps to categaorize and summarize genes by what they do. Two common ways to do this are GO terms and KEGG pathways.


```r
summary(res, alpha = 0.05)
```

```
## 
## out of 6041 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)     : 598, 9.9% 
## LFC < 0 (down)   : 1049, 17% 
## outliers [1]     : 6, 0.099% 
## low counts [2]   : 235, 3.9% 
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

## GO term enrichment

"A GO annotation is a statement about the function of a particular gene. GO annotations are created by associating a gene or gene product with a GO term. Together, these statements comprise a “snapshot” of current biological knowledge. Hence, GO annotations capture statements about how a gene functions at the molecular level, where in the cell it functions, and what biological processes (pathways, programs) it helps to carry out.

Different pieces of knowledge regarding gene function may be established to different degrees, which is why each GO annotation always refers to the evidence upon which it is based. All GO annotations are ultimately supported by the scientific literature, either directly or indirectly. In GO, the supporting evidence is presented in the form of a GO Evidence Codes and either a published reference or description of the methodology used to create the annotation. The GO evidence codes describe the type of evidence and reflect how far removed the annotated assertion is from direct experimental evidence, and whether this evidence was reviewed by an expert biocurator."  -- http://geneontology.org/docs/go-annotations/



```r
GO_df = toTable(org.Sc.sgdGO)
head(GO_df)
```

```
##   systematic_name      go_id Evidence Ontology
## 1            AWA1 GO:0007155      ISS       BP
## 2            ENA6 GO:0055085      ISA       BP
## 3            ENA6 GO:0055085      IMP       BP
## 4            ENA6 GO:0006814      ISA       BP
## 5            ENA6 GO:0006814      IMP       BP
## 6            ENS2 GO:0006314      IMP       BP
```

This frame comes with all three types of GO terms in one frame, BP = Biological Process, MF = Molecular Function, CC = Cellular COmponenet  


Convert the df to the format required for GOstats:

```r
goframeData = data.frame(GO_df$go_id, GO_df$Evidence, GO_df$systematic_name)
names(goframeData) = c("GO", "Evidence", "gene_id")
goframeData$GO <- as.character(goframeData$GO)
goframeData$Evidence <- as.character(goframeData$Evidence)
goframeData$gene_id <- as.character((goframeData$gene_id))
head(goframeData)
```

```
##           GO Evidence gene_id
## 1 GO:0007155      ISS    AWA1
## 2 GO:0055085      ISA    ENA6
## 3 GO:0055085      IMP    ENA6
## 4 GO:0006814      ISA    ENA6
## 5 GO:0006814      IMP    ENA6
## 6 GO:0006314      IMP    ENS2
```

Now turn this into a GO frame for GOstats:




















