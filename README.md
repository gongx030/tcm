# tcm

TCM is an R package for visualizing temporal single cell RNA-seq data.

## 1. Installation

The recommended installation method for `tcm` is using `install_github` command from `devtools` library.  You will first have to have [devtools](https://github.com/hadley/devtools) package installed.

```r
library(devtools)
install_github('gongx030/tcm')
```

A number of needed packages are installed in this process.

## 2. Quick Start

We first load the tcm package:
```r
library(tcm)
```

### 2.1 Visualizing a simulated temporal scRNA-seq dataset with 2,000 genes, 500 cells and five lineages.

Let us first simulate a simple temporal scRNA-seq data with 2,000 genes, 500 cells and five different lineages.  The single cell data are sampled across five time points following a sequentail differentiation model. The dropout noise was added using an exponential decay model. 

```r
set.seed(144)
sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, n.time.points = 5)
X <- assays(sim)$count 
time.table <- colData(sim)$time.table 
```

Note that X is the simulated gene by cell read count matrix (2,000 by 500), and time.table is a cell by time point binary matrix (500 by 5), indicating the timestamp of each cell.  

We run TCM on the simulated temporal scRNA-seq. 
```r
mf <- tcm(X, time.table = time.table)
```

We then plot the dimension reduction results from TCM:
```r
col.lineage <- rainbow(5) # color of cells from each lineage
bg.cell <- col.lineage[colData(sim)$lineage] # color of each cell
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
legend(par('usr')[2], par('usr')[4], 1:5, bty = 'n', xpd = NA, pt.bg = col.lineage, pch = 21, col = 'black', cex = 1.75)
```

![alt text](/docs/images/tcm_sim.png)

Let us the performance of t-SNE and diffusion map on the same simulated dataset: 

We first scale and log-transform the raw read counts data.
```r
X.log <- scale(log(X + 1))
```

We visualize the simulated temporal scRNA-seq data by t-SNE:
```r
library(Rtsne)
y <- Rtsne(t(X.log))$Y
plot(y[, 1], y[, 2], pch = 21, cex = 1.5, bg = bg.cell, col = 'black', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = 'tsne')
```
![alt text](/docs/images/tsne_sim.png)


Similarly, we visualize the dimension reduction results from diffusion map:
```r
library(diffusionMap)
y <- diffuse(dist(t(X.log)))$X[, c(1, 2)]
plot(y[, 1], y[, 2], pch = 21, cex = 1.5, bg = bg.cell, col = 'black', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = 'diffusion map')
```
![alt text](/docs/images/dm_sim.png)

We find that for a relatively easy dataset, both t-SNE and diffusion map are able to separate three lineages.  However, the visualization produced by TCM show much improved lineage trajectories. 

### 2.2 Visualizing a large-scale simulated temporal scRNA-seq dataset with 2,000 genes, 10,000 cells and five lineages.

```r
set.seed(2)
sim <- sim.rnaseq.ts(N = 2000, M = 10000, n.lineage = 5, n.time.points = 3)
X <- assays(sim)$count
time.table <- colData(sim)$time.table
```

We will use stochastic variational inference for optimizaing TCM, with batch size of 2,000 cells and using four CPU cores. It will take ~15 mins on an Intel Xeon 2.6GHz computer. 

```r
mf <- tcm(X, time.table = time.table, control = list(optimization.method = 'stochastic', batch.size = 2000, mc.cores = 4))
```

Visualization of the dimension reduction results:
```r
col.lineage <- rainbow(5)
bg.cell <- col.lineage[colData(sim)$lineage]
plot(mf, pch = 21, bg = bg.cell, cex = 0.5)
legend(par('usr')[2], par('usr')[4], 1:5, bty = 'n', xpd = NA, pt.bg = col.lineage, pch = 21, col = 'black', cex = 1.75)
```
![alt text](/docs/images/tcm_sim_large.png)


# 3. Case studies

## 3.1 Guo(2015), single cell PCR of early mouse embryonic development

```r
library(scDatasets)
data(guo)
X <- assays(guo)$count
time.table <- colData(guo)$time.table
X <- preprocess(X)

set.seed(1) 
mf <- tcm(X, time.table = time.table)
cg <- factor(colnames(time.table)[max.col(time.table)], colnames(time.table))
col.cg <- rainbow(nlevels(cg))
names(col.cg) <- levels(cg)
bg.cell <- col.cg[as.numeric(cg)]
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
```
![alt text](/docs/images/tcm_guo.png)

## 3.2 Single cell RNA-seq of the development of human primordial germ cells (PGC) and neighboring somatic cells from weeks 4 to 19 post-gestation. 
```r
library(scDatasets)
data(guo2)
X <- assays(guo2)$count
X <- preprocess(X)
time.table <- colData(guo2)$time.table

cg <- factor(colData(guo2)[['source.name']])
col.cg <- rainbow(nlevels(cg))
names(col.cg) <- levels(cg)
bg.cell <- col.cg[as.numeric(cg)]
bg.cell[colData(guo2)[['developmental.stage']] == '17 week gestation' & colData(guo2)[['gender']] == 'female' & colData(guo2)[['source.name']] == 'Primordial Germ Cells'] <- 'gold'
bg.cell[colData(guo2)[['developmental.stage']] == '19 week gestation' & colData(guo2)[['gender']] == 'male' & colData(guo2)[['source.name']] == 'Primordial Germ Cells'] <- 'purple'

set.seed(1)
mf <- tcm(X, time.table = time.table)
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
```
![alt text](/docs/images/tcm_guo2.png)



## 3.2 Single cell RNA-seq of iPSC-CM differentiation
```r
library(scDatasets)
data(iPSC_CM)
X <- preprocess2(assay(iPSC_CM), max.expressed.ratio = 1, min.expressed.gene = 0, min.expressed.cell = 2, normalize.by.size.effect = FALSE)
set.seed(1)
time.table <- colData(guo)[['time.table']]
mf <- tcm(X, K = 15, time = time.table, n.circle = 10, n.metacell = 15, n.prev = 2, remove.first.pc = TRUE, max.iter = 200)

cg <- factor(colnames(mf$V), c('D0', 'D6', 'D10', 'D30', 'D60'))
col.cg <- rainbow(nlevels(cg))
names(col.cg) <- levels(cg)
bg.cell <- col.cg[as.numeric(cg)]

dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
legend(par('usr')[2], par('usr')[4], levels(cg), bty = 'n', xpd = NA, pt.bg = col.cg, pch = 21, col = 'black')
```
![alt text](/docs/images/tcm_iPSC.png)

### Investigate the expression pattern of marker genes during iPSC-CM differentiation:

TNNT2 (cardiomyocyte marker)
```r
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = num2color(log(X['TNNT2', ] + 1)), cex = 2.25, main = 'TNNT2')
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/iPSC_TNNT2.png)

CDH11 (cardiofibroblast marker)
```r
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = num2color(log(X['CDH11', ] + 1)), cex = 2.25, main = 'CDH11')
```
![alt text](/docs/images/iPSC_CDH11.png)

## 3.3 Single cell RNA-seq of hESC derived mesodermal lineages
```r
library(scDatasets)
data(loh)
X <- preprocess2(assay(loh), max.expressed.ratio = 1, min.expressed.gene = 0, min.expressed.cell = 2, normalize.by.size.effect = FALSE)
set.seed(1)
time.table <- colData(loh)[['time.table']]
mf <- tcm(X, K = 10, time = time.table, n.circle = 10, n.metacell = 15, n.prev = 3, remove.first.pc = TRUE, max.iter = 50)

cn <- colData(loh)[['sample_alias']]  # cell names
lab2name <- c(
  'H7hESC' = 'hESC(D0)',
  'H7_derived_APS' = 'Anterior PS(D1)',
  'H7_derived_DLL1pPXM' = 'Paraxial Mesoderm(D2)',
  'H7_dreived_D2.25_Smtmrs' = 'Somitomeres(D2.25)',
  'H7_derived_ESMT' = 'Early Somite(D3)',
  'H7_derived_MPS' = 'Mid PS(D1)',
  'H7_derived_D2LtM' = 'Lateral Mesoderm(D2)',
  'H7_derived_D3GARPpCrdcM' = 'Cardiac Mesoderm(D3-D4)'
)
cg <- factor(cn, names(lab2name), lab2name)
col.cg <- rainbow(nlevels(cg))
names(col.cg) <- levels(cg)
bg.cell <- col.cg[as.numeric(cg)]
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
legend(par('usr')[2], par('usr')[4], levels(cg), bty = 'n', xpd = NA, pt.bg = col.cg, pch = 21, col = 'black')
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/tcm_Loh.png)

## 3.4 Single cell RNA-seq of human myoblast differentiation
```r
library(scDatasets)
data(trapnell)
X <- preprocess2(assay(trapnell), max.expressed.ratio = 1, min.expressed.gene = 0, min.expressed.cell = 2, normalize.by.size.effect = FALSE)
set.seed(1)
time.table <- colData(trapnell)[['time.table']]
mf <- tcm(X, K = 10, time = time.table, n.circle = 10, n.metacell = 15, n.prev = 3, remove.first.pc = TRUE, max.iter = 200)

cg <- colData(trapnell)[['time']]
col.cg <- rainbow(nlevels(cg))
names(col.cg) <- levels(cg)
bg.cell <- col.cg[as.numeric(cg)]
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
legend(par('usr')[2], par('usr')[4], levels(cg), bty = 'n', xpd = NA, pt.bg = col.cg, pch = 21, col = 'black')
```

TNNT2 (cardiac marker)
```r
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = num2color(log(X['ENSG00000118194', ] + 1)), cex = 2.25, main = 'TNNT2')
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/Trapnell_TNNT2.png)

PDGFRA (MC marker)
```r
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = num2color(log(X['ENSG00000134853', ] + 1)), cex = 2.25, main = 'PDGFRA')
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/Trapnell_PDGFRA.png)

IGFBP5 (skeletal marker)
```r
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = num2color(log(X['ENSG00000115461', ] + 1)), cex = 2.25, main = 'IGFBP5')
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/Trapnell_IGFBP5.png)


# 4. Session Information
```r
> sessionInfo()
R version 3.4.3 (2017-11-30)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.9 (Final)

Matrix products: default
BLAS: /panfs/roc/msisoft/R/3.4.3/lib64/R/lib/libRblas.so
LAPACK: /panfs/roc/msisoft/R/3.4.3/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices datasets  utils     methods   base

other attached packages:
 [1] scDatasets_0.0.3           wordcloud_2.5              RColorBrewer_1.1-2         plotrix_3.7                gplots_3.0.1               igraph_1.2.1               MASS_7.3-47                cluster_2.0.6              FNN_1.1                    gtools_3.5.0               fields_9.6                 maps_3.2.0
[13] spam_2.1-2                 dotCall64_0.9-5.2          irlba_2.3.2                Matrix_1.2-12              SummarizedExperiment_1.8.1 DelayedArray_0.4.1         matrixStats_0.53.1         Biobase_2.38.0             GenomicRanges_1.30.3       GenomeInfoDb_1.14.0        IRanges_2.12.0             S4Vectors_0.16.0
[25] BiocGenerics_0.24.0        BiocInstaller_1.28.0

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.16           compiler_3.4.3         XVector_0.18.0         bitops_1.0-6           tools_3.4.3            zlibbioc_1.24.0        digest_0.6.15          memoise_1.1.0          lattice_0.20-35        pkgconfig_2.0.1        commonmark_1.4         GenomeInfoDbData_1.0.0 withr_2.1.2            stringr_1.3.0
[15] roxygen2_6.0.1         xml2_1.2.0             caTools_1.17.1         devtools_1.13.5        R6_2.2.2               gdata_2.18.0           magrittr_1.5           KernSmooth_2.23-15     stringi_1.1.7          RCurl_1.95-4.10        slam_0.1-42
```
