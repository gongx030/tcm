% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcm.r
\name{plot.tcm}
\alias{plot.tcm}
\title{plot.tcm}
\usage{
\method{plot}{tcm}(x, ...)
}
\arguments{
\item{x}{a tcm object}

\item{...}{Further graphical parameters may also be supplied as arguments}
}
\description{
Plot the TCM results
}
\examples{
library(tcm)
set.seed(122)
sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, n.time.points = 5)
X <- assays(sim)$count
time.table <- colData(sim)$time.table
mf <- tcm(X, time.table = time.table)
bg.cell <- rainbow(5)[colData(sim)$lineage]
dev.new(height = 10, width = 12)
par(mar = c(5, 5, 5, 15))
plot(mf, pch = 21, bg = bg.cell, cex = 1.5)

}
\author{
Wuming Gong, \email{gongx030@umn.edu}
}
