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

### 2.1 A simulated temporal scRNA-seq dataset with 2,000 genes, 500 cells and five lineages.

Let us first simulate a simple temporal scRNA-seq data with 2,000 genes, 500 cells and five different lineages.  The single cell data are sampled across five time points following a sequentail differentiation model. The dropout noise was added using an exponential decay model. 

```r
set.seed(16)
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
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
legend(par('usr')[2], par('usr')[4], 1:5, bty = 'n', xpd = NA, pt.bg = col.lineage, pch = 21, col = 'black', cex = 1.75)```

We can visualize the inferred developmental trajectory:
```r
mf <- trajectory(mf)
add.paths(mf, lwd = 2)
```
![alt text](/docs/images/tcm_sim.png)

Let us see whether or not popular tools such as t-SNE and diffusion map are able separate three lineages: 

We first scale and log-transform the raw read counts data.
```r
X2 <- scale(log(sim$X + 1))
```

We visualize the simulated temporal scRNA-seq data by t-SNE:
```r
library(Rtsne)
set.seed(1)
y.tsne <- t(Rtsne(t(X2), check_duplicates = FALSE)$Y)
dev.new(height = 8, width = 10)
par(mar = c(5, 5, 5, 15))
plot(y.tsne[1, ], y.tsne[2, ], pch = 21, cex = 1.75, bg = bg.cell, col = 'black', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = 't-SNE')
legend(par('usr')[2], par('usr')[4], 1:sim$n.lineage, bty = 'n', xpd = NA, pt.bg = bg.lineage, pch = 21, col = 'black', cex = 1.75)
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/tsne_sim0.png)


Similarly, we visualize the dimension reduction results from diffusion map:
```r
library(diffusionMap)
set.seed(1)
y <- t(diffuse(dist(t(X2)))$X)[1:2, ]
dev.new(height = 8, width = 10)
par(mar = c(5, 5, 5, 15))
plot(y[1, ], y[2, ], pch = 21, cex = 1.5, bg = bg.cell, col = 'black', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = 'diffusion map')
legend(par('usr')[2], par('usr')[4], 1:sim$n.lineage, bty = 'n', xpd = NA, pt.bg = bg.lineage, pch = 21, col = 'black', cex = 1.75)
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/dm_sim0.png)

We find that for a relatively easy dataset, both t-SNE and diffusion map are able to separate three lineages.  However, the visualization produced by TCM show much improved lineage trajectories. 

### 2.2 A simulated temporal scRNA-seq dataset with 2,000 genes, 500 cells and five lineages.

Let us then turn to a more heterogenous dataset.  We simulate a temporal scRNA-seq dataset with 2,000 genes, 500 cells and **five different lineages**. The single cell data are sampled across **five time points** following a sequential differentiation model and an exponential decay model for the dropout noise. 

```r
set.seed(1)
sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, type = 'sequential', n.time.points = 5)
```

We visualize dimension reduction results produced by TCM:
```r
time.table <- table(1:sim$M, factor(sim$time))
set.seed(1)
mf <- tcm(sim$X, K = 15, time = time.table, n.circle = 10, n.metacell = 15, n.prev = 3, max.iter = 50)
bg.lineage <- rainbow(sim$n.lineage) # color for each lineage
bg.cell <- bg.lineage[sim$lineage] # color for each cell
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
legend(par('usr')[2], par('usr')[4], 1:sim$n.lineage, bty = 'n', xpd = NA, pt.bg = bg.lineage, pch = 21, col = 'black', cex = 1.75)
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/tcm_sim.png)


We scale and log-transform the raw read counts data. 
```r
X2 <- scale(log(sim$X + 1))
```

Let us visualize the simulated temporal scRNA-seq data by t-SNE and diffusion map:
```r
set.seed(1)
y.tsne <- t(Rtsne(t(X2), check_duplicates = FALSE)$Y)
dev.new(height = 8, width = 10)
par(mar = c(5, 5, 5, 15))
plot(y.tsne[1, ], y.tsne[2, ], pch = 21, cex = 1.75, bg = bg.cell, col = 'black', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = 't-SNE')
legend(par('usr')[2], par('usr')[4], 1:sim$n.lineage, bty = 'n', xpd = NA, pt.bg = bg.lineage, pch = 21, col = 'black', cex = 1.75)
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/tsne_sim.png)

```r
set.seed(1)
y <- t(diffuse(dist(t(X2)))$X)[1:2, ]
dev.new(height = 8, width = 10)
par(mar = c(5, 5, 5, 15))
plot(y[1, ], y[2, ], pch = 21, cex = 1.5, bg = bg.cell, col = 'black', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = 'diffusion map')
legend(par('usr')[2], par('usr')[4], 1:sim$n.lineage, bty = 'n', xpd = NA, pt.bg = bg.lineage, pch = 21, col = 'black', cex = 1.75)
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/dm_sim.png)

This example clearly demonstrates TCM has significantly better capability of revealing the hidden lineages from very heterogenous temporal scRNA-seq datasets than traditional dimension reduction methods such as t-SNE and diffusion map. 

# 3. Case studies

## 3.1 Guo(2015), single cell PCR of early mouse embryonic development

```r
library(scDatasets)
data(guo)
X <- preprocess2(assay(guo), max.expressed.ratio = 1, min.expressed.gene = 0, min.expressed.cell = 2)
set.seed(1)
time.table <- colData(guo)[['time.table']]
mf <- tcm(X, K = 10, time = time.table, n.circle = 10, n.metacell = 15, n.prev = 3, max.iter = 50)
cg <- factor(colnames(time.table)[max.col(time.table)], colnames(time.table))
col.cg <- rainbow(nlevels(cg))
names(col.cg) <- levels(cg)
bg.cell <- col.cg[as.numeric(cg)]
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/tcm_Guo.png)

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
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/tcm_iPSC.png)

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
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/iPSC_CDH11.png)

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

## 3.5 Single cell RNA-seq of the development of human primordial germ cells (PGC) and neighboring somatic cells from weeks 4 to 19 post-gestation. 
```r
library(scDatasets)
data(guo2)
X <- preprocess2(assay(guo2), max.expressed.ratio = 1, min.expressed.gene = 0, min.expressed.cell = 2, normalize.by.size.effect = FALSE)
set.seed(1)
time.table <- colData(guo2)[['time.table']]
mf <- tcm(X, K = 15, time = time.table, n.circle = 10, n.metacell = 15, n.prev = 3, remove.first.pc = TRUE, max.iter = 200)

cg <- factor(colData(guo2)[['source.name']])
col.cg <- rainbow(nlevels(cg))
names(col.cg) <- levels(cg)
bg.cell <- col.cg[as.numeric(cg)]
bg.cell[colData(guo2)[['developmental.stage']] == '17 week gestation' & colData(guo2)[['gender']] == 'female' & colData(guo2)[['source.name']] == 'Primordial Germ Cells'] <- 'gold'
bg.cell[colData(guo2)[['developmental.stage']] == '19 week gestation' & colData(guo2)[['gender']] == 'male' & colData(guo2)[['source.name']] == 'Primordial Germ Cells'] <- 'purple'
ab2pch <- c('male' = 21, 'female' = 22)
pch.cell <- lab2pch[colData(guo2)[['gender']]]

dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
legend(par('usr')[2], par('usr')[4], levels(cg), bty = 'n', xpd = NA, pt.bg = col.cg, pch = 21, col = 'black')
```
![alt text](https://github.com/gongx030/tcm/blob/master/docs/images/tcm_Guo2.png)


# 4. Session Information
```r
> sessionInfo()
R version 3.3.3 (2017-03-06)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.9 (Final)

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base

other attached packages:
[1] tcm_1.0.2            BiocInstaller_1.24.0

loaded via a namespace (and not attached):
 [1] lattice_0.20-35    FNN_1.1            gtools_3.5.0       bitops_1.0-6       MASS_7.3-47        grid_3.3.3         magrittr_1.5       spam_2.1-1         KernSmooth_2.23-15 gplots_3.0.1       irlba_2.3.1        gdata_2.18.0       Matrix_1.2-11      igraph_1.1.2       maps_3.2.0         fields_9.0
[17] plotrix_3.6-6      pkgconfig_2.0.1    cluster_2.0.6      caTools_1.17.1     dotCall64_0.9-04
```
