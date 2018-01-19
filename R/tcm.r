library(irlba)	# irlba
library(fields)	# rdist
library(gtools)	# ddirichlet
library(FNN)	# knnx.index
library(cluster) # pam
library(MASS) # mvrnorm
library(Matrix)
library(igraph)	# get.shortest.paths
library(gplots)	# colorpanel
library(plotrix)	# draw.circle
library(SummarizedExperiment) 	
library(parallel) # mclapply


#' Simulating time series single cell RNA-seq data
#' 
#' @param N number of genes
#' @param M number of cells
#' @param K number of metagenes
#' @param n.circle number of circles on the latent space
#' @param n.circle prototypes per circle on the latent space
#' @param lambda the argument of a exponential decay model for introducing the dropout noise
#' @param alpha0 the argument for a Dirichlet distribution for sampling metagene basis
#' @param library.size total number of reads per cell
#' @param n.neighbors number of neighbors when looking for nearby prototypes
#' @param n.lineage number of simulated lineages
#' @param type the type of differentiation models
#' @param n.time.points the number of simulated time points
#' 
#' @return A SummarizedExperiment object of simulated temporal scRNA-seq data
#'
#' @author Wuming Gong, \email{gongx030@umn.edu}
#' 
#' @examples
#' set.seed(1)
#' sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, type = 'sequential', n.time.points = 5)
#' 
#' @export
#'
#' @importFrom fields rdist
#' @importFrom irlba irlba
#' @importFrom gtools ddirichlet rdirichlet
#' @importFrom cluster pam
#' @importFrom FNN knnx.index
#' @importFrom MASS mvrnorm
#' @import Matrix
#' @importFrom graphics plot points
#' @importFrom stats as.dist cmdscale rbinom rmultinom rnorm runif
#' @importFrom igraph get.shortest.paths graph.adjacency
#' @importFrom gplots colorpanel
#' @importFrom plotrix draw.circle
#' 
sim.rnaseq.ts <- function(N = 200, M = 50, K = 20, n.circle = 100, n.metacell = 10, lambda = 0.25, alpha0 = 0.1, library.size = c(1e3, 1e5), n.neighbors = 7, n.lineage = 3, choose.leaves = 'random', type = 'sequential', n.time.points = 10){

	# generate a metacell landscape
	C <- expand.grid(angle = seq(0, 2 * pi, length.out = n.metacell + 1)[-1], r = seq(0, 1, length.out = n.circle + 1)[-1])	# polar coord
	C <- cbind(C, x = C[, 'r'] * cos(C[, 'angle']), y = C[, 'r'] * sin(C[, 'angle']))	# x-y coord
	C <- transform(C, circle = as.numeric(factor(r)), pos = as.numeric(factor(angle)))
	C <- cbind(C, data.frame(active = FALSE, parent = NA))
	C <- rbind(data.frame(angle = 0, r = 0, x = 0, y = 0, circle = 0, pos = 0, active = FALSE, parent = NA), C)
	C <- cbind(C, index = 1:nrow(C))

	H0 <- nrow(C)	# number of metacells

	L.inv <- rdist(as.matrix(C[, c('x', 'y')]))
	L.inv <- exp(-L.inv^2 / 2)	# squared exponential covariance function
	Theta <- mvrnorm(K, rep(0, H0), Sigma = L.inv)	# randomly generate the metacell coefficients

	A <- NULL
	for (i in 1:n.circle){
		Cc <- C[C[, 'circle'] == i, , drop = FALSE]
		Cp <- C[C[, 'circle'] == i - 1, , drop = FALSE]
		nn <- knnx.index(Cc[, c('x', 'y')], Cp[, c('x', 'y')], n.neighbors)
		A <- rbind(A, data.frame(from = Cp[rep(1:nrow(nn), n.neighbors), 'index'], to = Cc[c(nn), 'index']))
	}
	A <- sparseMatrix(i = A[, 1], j = A[, 2], x = sqrt(colSums((Theta[, A[, 1]] - Theta[, A[, 2]])^2)), dims = c(H0, H0))
	A <- (A + t(A)) / 2
	A <- graph.adjacency(A, weighted = TRUE, diag = FALSE)

	if (choose.leaves == 'random')
		leaves <- sample(which(C[, 'circle'] == n.circle), n.lineage)
	else if (choose.leaves == 'equal'){
		hs <- which(C[, 'circle'] == n.circle)
		leaves <- round(seq(min(hs), max(hs), length.out = n.lineage + 1))[-1]
	}

	paths <- lapply(leaves, function(leaf) get.shortest.paths(A, from = 1, to = leaf)$vpath[[1]])
	C[unlist(paths), 'active'] <- TRUE
	LN <- matrix(FALSE, H0, n.lineage)	# active metacells for each lineage
	for (i in 1:length(paths)) LN[paths[[i]], i] <- TRUE

	H <- sum(C[, 'active'])	# number of active cellular states

	if (type == 'sequential'){
		starts <- round(seq(1, n.circle, length.out = n.time.points + 1))[-(n.time.points + 1)]	# start of time interval
		ends <- c(starts[-1] - 1, n.circle)
	}else if (type == 'delayed'){
		starts <- rep(1, n.time.points)
		ends <- round(seq(1, n.circle, length.out = n.time.points + 1))[-1]
	}else if (type == 'forward'){
		starts <- round(seq(1, n.circle, length.out = n.time.points + 1))[-(n.time.points + 1)]	# start of time interval
		ends <- rep(n.circle, n.time.points)
	}
	P <- do.call('cbind', lapply(1:n.time.points, function(i) 1:n.circle >= starts[i] & 1:n.circle <= ends[i]))	# n.circle ~ n.time.points
	P <- P %*% Diagonal(x = 1 / Matrix::colSums(P))

	GB <- sparseMatrix(i = 2:H0, j = C[-1, 'circle'], dims = c(H0, n.circle))	# metacell ~ layers, excluding the root
	GBP <- GB %*% P	# metacell ~ time, sampling probability

	z <- rep(NA, M)	# metacell assignment for each cell
	lineage <- rep(NA, M)	# lineage assignment for each cell
	time <- rep(NA, M)	# time assignment for each cell
	for (m in 1:M){
		lineage[m] <- sample(1:n.lineage, 1)	# randomly assign a lineage
		time[m] <- sample(1:n.time.points, 1)	# randomly assign a time point
		z[m] <- sample(which(LN[, lineage[m]]), 1, prob = GBP[LN[, lineage[m]], time[m]])
	}

	V <- Theta[, z]	# metacell coefficient for each cell

	V.exp <- exp(V) %*% diag(1 / colSums(exp(V)))
	libsize <- runif(M, library.size[1], library.size[2])	# sampling the library size for each cell
	U <- t(rdirichlet(K, alpha = rep(alpha0, N)))	# sampling the metagene basis
	Mu <- do.call('cbind', lapply(1:M, function(m) rmultinom(1, libsize[m], prob = U %*% V.exp[, m, drop = FALSE])))	# sampling the number of reads
	prob <- exp(-lambda * log(Mu + 1))	# introduce the dropout noise by a exponential decay model
	D <- matrix(rbinom(N * M, 1, prob), nrow = N, ncol = M) == 1
	X <- Mu
	X[D] <- 0
	rownames(X) <- rownames(Mu) <- sprintf('G%d', 1:N)

	SummarizedExperiment(
		assays = list(count = X, Mu = as.matrix(Mu), D = D), 
		rowData = DataFrame(U = I(U)),
		colData = DataFrame(V = I(t(V)), z = z, lineage = lineage, circle = C[z, 'circle'], time = I(time), time.table = I(table(1:M, factor(time)))),
		metadata = list(N = N, M = M, Theta = Theta, C = C, paths = paths, n.circle = n.circle, n.metacell = n.metacell, A = A, n.lineage = n.lineage)
	)

} # end of sim.rnaseq.ts


#' Visualizing temporal scRNA-seq by topographic cell map 
#' 
#' @param X a normalized read count matrix where each row represents a gene and each column represents a cell
#' @param K the number of metagenes
#' @param time.table a cell by time point table indicating the source of each cell
#' @param landscape the prototype landscape for TCM. The landscape need to be initialized by init.landscape(). 
#' @param optimization.method the method for optimizing TCM, either 'batch' or 'stochastic' (default: stochastic)
#' @param batch.size the size of randomly sampled cells for updating global variables in stochastic variational inference (default: 100).
#'				Note that batch.size is ignored if 'optimization.method' is 'batch'
#' @param max.iter maximum iterations (default: 100)
#' @param mc.cores # of CPU cores for optimization (default: 1)
#' 
#' @return a tcm object
#
#' @author Wuming Gong, \email{gongx030@umn.edu}
#
#' @examples
#' set.seed(1)
#' sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, type = 'sequential', n.time.points = 5)
#' bg.lineage <- rainbow(sim$n.lineage)
#' bg.cell <- bg.lineage[sim$lineage]
#' time.table <- table(1:sim$M, factor(sim$time))
#' set.seed(1)
#' mf <- tcm(sim$X, time = time.table, n.circle = 10, n.metacell = 15, n.prev = 3, max.iter = 10)
#' set.seed(1)
#' dev.new(height = 10, width = 12)
#' par(mar = c(5, 5, 5, 15))
#' plot(mf, pch = 21, bg = bg.cell, cex = 2.25)
#' legend(par('usr')[2], par('usr')[4], 1:sim$n.lineage, bty = 'n', xpd = NA, pt.bg = bg.lineage, pch = 21, col = 'black', cex = 1.75)
#'
#' @export
#' 
tcm <- function(X, time.table = NULL, landscape = NULL, nc = 3, init.V.method = 'mds', init.U.method = 'pam', optimization.method = 'stochastic', batch.size = 2000, max.iter = 100, mc.cores = 1){

	N <- nrow(X)
	M <- ncol(X)
	eps <- .Machine$double.eps

	if (is.null(time.table) && landscape[['type']] == 'temporal.convolving')
		stop('time.table cannot be NULL if landscape$type is temporal.convolving')
	
	if (is.null(landscape))
		stop('landscape cannot be NULL')

	K <- landscape$K

	CT <- as(as(time.table, 'matrix'), 'ngCMatrix') # cell ~ time
	time <- colnames(time.table)

	cat('-------------------------------------------------------------------------------------------------\n')
	cat('topographic cell map\n')
	cat('-------------------------------------------------------------------------------------------------\n')
	cat(sprintf('[%s] number of input genes(N): %d\n', Sys.time(), N))
	cat(sprintf('[%s] number of input cells(M): %d\n', Sys.time(), M))
	cat(sprintf('[%s] optimization method: %s\n', Sys.time(), optimization.method))
	cat(sprintf('[%s] method for initializing U: %s\n', Sys.time(), init.U.method))
	cat(sprintf('[%s] method for initializing V: %s\n', Sys.time(), init.V.method))
	if (optimization.method == 'stochastic'){
		cat(sprintf('[%s] batch.size: %d\n', Sys.time(), batch.size))
	}
	cat(sprintf('[%s] number of cores: %d\n', Sys.time(), mc.cores))

#	gs <- max.col(time.table)
#	pvalues <- unlist(mclapply(1:N, function(n) kruskal.test(split(X[n, ], list(gs)))$p.value, mc.cores = mc.cores))
#	fdr <- p.adjust(pvalues, method = 'fdr')
#	X <- X[fdr < 0.001, ]
#	N <- nrow(X)
#	cat(sprintf('[%s] number of temporally dynamically expressed genes (adjusted kruskal test p<0.05): %d\n', Sys.time(), N))

	cat(sprintf('[%s] initializing metagene coefficients and basis\n', Sys.time()))

	if (init.V.method == 'irlba'){
		V <- scale(t(irlba(scale(log(X + 1)), nu = 1, nv = K + 1)$v[, 2:(K + 1)]))
	}else if (init.V.method == 'mds'){
		X.log <- scale(log(X + 1))	# scaled and log transformed input
		s <- irlba(X.log, nu = 1, nv = 1)
		X.log <- scale(X.log - (s$u[, 1] %o% s$v[, 1]) * s$d)
		D <- rdist(t(X.log))^2	# cell-wise distance matrix
		V <- scale(t(cmdscale(as.dist(D), eig = TRUE, k = K)$points))
	}
	V <- matrix(c(V), nrow(V), ncol(V))

	a <- rep(0, M)	# initialize the cell effect
	if (init.U.method == 'random'){
		U <- matrix(runif(N * K), N, K) 	# initialize the metagene basis
	}else if (init.U.method == 'pam'){
		CC <- sparseMatrix(i = 1:M, j = pam(t(V), K, cluster.only = TRUE), dims = c(M, K))	# # cells ~ clusters
		U <- as.matrix(X %*% CC %*% diag(1 / Matrix::colSums(CC)))	# initialize E[U]
	}

	mf <- ccr(X, U = U, V = V, a = a, optimization.method = optimization.method, batch.size = batch.size, max.iter = max.iter, mc.cores = mc.cores)
	mf <- ccm(X, U = mf$U, V = V, a = mf$a, optimization.method = optimization.method, batch.size = batch.size, max.iter = max.iter, mc.cores = mc.cores)
	V <- mf$V	# the metagene coefficient
	U <- mf$U	# the metagene basis
	a <- mf$a

	cat(sprintf('[%s] initializing prototype landscape\n', Sys.time()))
	if (landscape[['type']] == 'temporal.convolving'){	# initialization for landscape type = temporal.convolving

		ls <- init.landscape(type = 'plate', K = K, n.prototype = n.prototype, n.circle = 5)
		mf <- gtm(V = V, landscape = ls, optimization.method = optimization.method, min.beta = 1, batch.size = batch.size, max.iter = max.iter, mc.cores = mc.cores)
		CC <- sparseMatrix(i = 1:M, j = max.col(mf$Z), dims = c(M, ls[['H.prototype']])) %*% ls[['MC']]	# cell ~ circle
		CC <- CC[, Matrix::colSums(CC) > 0]	# circles with mapped cells

		for (i in ncol(CC):1){	# from outer to inner region
			i <- ncol(CC)
			m <- which(Matrix::rowSums(CC[, i:ncol(CC), drop = FALSE]) > 0)	# cells mapped to current circle
			mf <- gtm(V = V[, m], landscape = landscape, CT = CT[m, ], optimization.method = optimization.method, min.beta = 0, batch.size = batch.size, max.iter = max.iter, mc.cores = mc.cores)
			landscape <- mf$landscape
		}

		mf <- tcm.core(X = X, U = U, V = V, a = a, landscape = landscape, CT = CT, optimization.method = optimization.method, batch.size = batch.size, max.iter = max.iter, mc.cores = mc.cores)
	}

	mf

} # end of tcm


#' Computing y[i, j] = exp(x[i, j]) / sum(exp(x[, j]))
#'
softmax <- function(x){
	x.max <- apply(x, 2, max)
	y <- log(Matrix::rowSums(exp(t(x) - x.max)) + .Machine$double.eps) + x.max	# log(sum(exp(x[, j])))
	y <- t(exp(t(x) - y))
} # end of softmax 


#' The core function for optimizaing TCM
#'
#' @param	S	a column stochastic matrix (colSums(S) == 1).  describing the relationship between free and convolving prototypes
#' @param	L Precision matrix defined for the prototypes
#' @param optimization.method optimization method for TCM, either 'batch' or 'stochastic'
#' @param batch.size size of mini batch
#' @param mc.cores number of CPU cores (default: 1)
#' @param max.size.per.batch maximum samples per batch
#'
tcm.core <- function(X, U = NULL, V = NULL, a = NULL, landscape = NULL, CT = NULL, optimization.method = 'batch', min.beta = 0, c0 = 0.001, d0 = 0.001, max.iter = 100, verbose = TRUE, batch.size = 100, decay.rate = 1, forgetting.rate = 0.75, mc.cores = 1, max.size.per.batch = 1000){

	eps <- .Machine$double.eps

	# the number of iterations for Newton step for optimizing V and V.var
	max.iter2 <- switch(optimization.method, 'batch' = 50, 'stochastic' = 1)

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- ncol(U)	# number of metagenes
	H <- landscape[['H.prototype']]	# number of prototypes

	if (!is.null(landscape)){
		Theta <- landscape[['Theta.free']] %*% landscape[['S']]	# all prototypes
		if (!is.null(CT))
			W <- CT %*% landscape[['TMC']]	# a binary cell ~ prototype assignment matrix
		else
			W <- matrix(TRUE, M, landscape[['H.prototype']])
	}else
		stop('landscape must be specified')

	if (is.null(a))
		a <- rep(0, M)	# the cell effect

	# splitting the samples into batches
	if (M <= mc.cores * max.size.per.batch){	
		n.batch <- mc.cores
	}else{
		n.batch <- ceiling(M / max.size.per.batch)
	}
	groups <- sample(1:n.batch, M, replace = TRUE)
	size <- table(factor(groups, 1:n.batch))	# number of samples per core

	# initialize the local variables
	localvar <- mclapply(1:n.batch, function(b){
		m <- groups == b
		V.var <- matrix(1, nrow = K, ncol = size[b])	# Var[V]
		V.exp <- exp(V[, m] + 1/2 * V.var)	# E[exp(V)]
		list(a = a[m], V = V[, m], V.exp = V.exp, Z = matrix(1 / H, size[b], H))
	}, mc.cores = mc.cores)

	iter <- 1
	optval <- NULL
	while (iter <= max.iter){

		Theta.free.p <- landscape[['Theta.free']]

		# updating all local variables: E[V], E[exp(V)], a and E[Z]
		localvar <- mclapply(1:n.batch, function(b){
			m <- groups == b
			update.localvar(
				X = X[, m, drop = FALSE],
				U = U,
				V = localvar[[b]]$V, 
				V.exp = localvar[[b]]$V.exp,
				W = W[m, ],
				landscape = landscape,
				model = 'tcm', max.iter = max.iter2
			)
		}, mc.cores = mc.cores)

		rho <- switch(optimization.method,
			'stochastic' = learning.rate(iter, decay.rate, forgetting.rate),
			'batch' = 1
		)

		if (optimization.method == 'stochastic'){
			# SVI: sampling a subset of samples(cells) for updating global variables
			m <- lapply(1:n.batch, function(b) sample(1:size[b], min(max(1, round(batch.size / n.batch)), size[b])))
		}else{
			# batch optimization: use all the samples
			m <- lapply(1:n.batch, function(b) 1:size[b])
		}

		# updating free prototypes Theta
		SZS <- landscape[['S']] %*% Diagonal(x = landscape[['beta']] * Reduce('+', mclapply(1:n.batch, function(b) colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]), mc.cores = mc.cores))) %*% t(landscape[['S']])
		SZSL.inv <- chol2inv(chol(as.matrix(SZS + landscape[['L']])))
		landscape[['Theta.free']] <- (1 - rho) * landscape[['Theta.free']] + rho * landscape[['beta']] * Reduce('+', mclapply(1:n.batch, function(b){
			localvar[[b]]$V[, m[[b]], drop = FALSE] %*% localvar[[b]]$Z[m[[b]], , drop = FALSE] %*% t(landscape[['S']]) %*% SZSL.inv
		}, mc.cores = mc.cores))

#		U <- (1 - rho) * U + rho * Reduce('+', mclapply(1:n.batch, function(b){
#			c <- U * (X[, which(groups == b)[m[[b]]]] / (U %*% localvar[[b]]$V.exp[, m[[b]], drop = FALSE] + eps)) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE]) + c0
#			d <- matrix(1, N, length(m[[b]])) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE] %*% diag(exp(localvar[[b]]$a[m[[b]]]))) + d0
#			U.log <- digamma(c) - log(d + eps)	# E[log(U)]
#			exp(U.log)
#		}, mc.cores = mc.cores)) / n.batch

		# updating beta
		landscape[['beta']] <- (1 - rho) * landscape[['beta']] + rho * Reduce('+', mclapply(1:n.batch, function(b){
			(K * sum(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + 0.001) / (sum(localvar[[b]]$Z[m[[b]], , drop = FALSE] * rdist(t(localvar[[b]]$V[, m[[b]], drop = FALSE]), t(Theta))^2) + 0.001)
			}, mc.cores = mc.cores)) / n.batch
		landscape[['beta']] <- max(landscape[['beta']], min.beta)

		# updating kappa/pi.log
		kappa <- Reduce('+', mclapply(1:n.batch, function(b){
			colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + landscape[['alpha0']]
		}, mc.cores = mc.cores)) 
		landscape[['pi.log']] <- (1 - rho) * landscape[['pi.log']] + rho * (digamma(kappa) - digamma(sum(kappa))) # E[log(pi)]
		

		if (iter == 1 || iter %% 10 == 0){
			J <- norm(landscape[['Theta.free']] - Theta.free.p, '2')
			n.cluster <- length(unique(unlist(mclapply(1:n.batch, function(b) max.col(localvar[[b]]$Z), mc.cores = mc.cores))))
			cat(sprintf('[%s] tcm | iter=%5.d | M=%5.d | n.batch=%3.d | dJ=%7.3e | # non-empty prototypes=%4.d | beta=%7.3f\n', Sys.time(), iter, M, n.batch, J, n.cluster, landscape[['beta']]))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), n.cluster = n.cluster, beta = landscape[['beta']], J = J))
			if (J < 0.01)
				break
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	a <- rep(0, M)
	V <- matrix(0, K, M)
	Z <- matrix(0, M, H)
	for (b in 1:mc.cores){
		V[, groups == b] <- localvar[[b]]$V
		Z[groups == b, ] <- localvar[[b]]$Z
		a[groups == b] <- localvar[[b]]$a
	}

	structure(list(H = H, N = N, M = M, K = K, U = U, V = V, a = a, Z = Z, CT = CT, landscape = landscape, optval = optval), class = 'tcm')

} # end of tcm.core



#' plot tcm object.
#' 
#' @param x tcm object
#' @param ... Further graphical parameters may also be supplied as arguments
#'
#' @export
#' 
plot.tcm <- function(x, ...){

	param <- list(...)

	mf <- x

	if (mf$landscape$type %in% c('temporal.convolving', 'plate')){

		H <- mf$landscape$H	# number of prototypes
		K <- mf$landscape$K	# number of metagenes
		M <- ncol(mf$V)	# number of cells
		mem <- max.col(mf$Z)	# cluster membership
		coord <- mf$landscape$Y.prototype
		csize <- table(factor(mem, 1:H))	# size of cells that were assigned to each prototype
		ncls <- sum(csize > 0)	# number of prototypes that have at least one assigned cells
		centers <- which(csize > 0)	# index of non-empty prototype
		CT <- mf$CT
		Yp <- mf$landscape$Yp.prototype
		Y <- mf$landscape$Y.prototype

		n.prototype <- mf$landscape$n.prototype

		lim <- c(-max(Yp[, 'r']), max(Yp[, 'r']))	# limit of the 2D space
		plot(NA, xlim = lim * 1.3, ylim = lim * 1.3, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', main = param[['main']])

		rs <- unique(Yp[, 'r'])
		col.time <- colorpanel(length(rs), low = 'gray100', mid = 'gray70', high = 'gray40')
		for (i in length(rs):1){
			draw.circle(0, 0, rs[i], nv = 100, border = NA, col = col.time[i])
		}

		if (mf[['landscape']][['type']] == 'temporal.convolving'){
			for (t in 1:ncol(mf[['CT']]))
				draw.circle(0, 0, t, nv = 100, border = 'black', lty = 2, lwd = 2)
		}

		Y.cell2 <- t(jitter2(t(coord[, max.col(mf$Z)]), amount = 0.125 * max(lim)))	# jittered coordinate of cells
		points(Y.cell2[1, ], Y.cell2[2, ], ...)

	}

} # end of plot.tcm


#' Adding an Enclidean noise, x must be in a two dimensional matrix
#'
jitter2 <- function(x, amount = 0.1){
	phi <- runif(nrow(x), min = 0, max = 2 * pi)
	r <- runif(nrow(x), min = 0, max = amount)
	x[, 1] <- x[, 1] + r * cos(phi)
	x[, 2] <- x[, 2] + r * sin(phi)
	x
} # end of jitter2


#' Scale a continuous vector into color representation
#'
num2color <- function(x, cp = colorpanel(100, low = 'blue', mid = 'black', high = 'yellow')) cp[round((x - min(x)) / (max(x) - min(x)) * (length(cp) - 1)) + 1]


#' learning rate for stochastic variational inference
#'
#' @param t current iteration step
#' @param decay.rate down-weights for early iterations (default: 5)
#' @param forgetting.rate controls how quickly old information is forgotten (default: 0.75)
#'
learning.rate <- function(t, decay.rate = 5, forgetting.rate = 0.75) ifelse(t == 0, 1, (t + decay.rate)^(-forgetting.rate))


#' updating local variable E[V], E[exp(V)], a and E[Z]
#' @param X a N by M expression matrix
#' @param U a N by K metagene basis matrix
#' @param V the initial values of a K by M metagene coefficient matrix
update.localvar <- function(X, U, V = NULL, V.exp = NULL, W = NULL, landscape = NULL, model = NA, max.iter = 1){

	eps <- .Machine$double.eps

	if (!is.null(W))
		W <- as.matrix(W)
	
	if(model == 'tcm' || model == 'ccm' || model == 'ccr'){

		K <- ncol(U)
		N <- nrow(U)
		M <- ncol(X)

		if (is.null(V.exp))
			V.exp <- exp(V)

		Xp <- U %*% V.exp + eps
		SN <- as.matrix(V.exp * ( t(U) %*% ( X / Xp) ))
		a <- log(colSums(SN)) - log(colSums(Xp) + eps)	# updating the cell effect (a)

		if (model == 'tcm'){
			B <- matrix(a, nrow = K, ncol = M, byrow = TRUE) + log(t(U) %*% matrix(1, N, M) + eps) # sum of a[m] + log(U[n, k] * W[n, m]) from n=1 to n=N
			Theta <- landscape[['Theta.free']] %*% landscape[['S']]	# all prototypes
			H <- landscape[['H.prototype']]	# number of prototypes
			# construct a permutation matrix to compute \sum_{m=1}^M \sum_{h=1}^H ||\Theta_h - V_m||^2
			perm <- expand.grid(1:M, (M + 1):(M + H))
			P <- sparseMatrix(i = c(perm[, 1], perm[, 2]), j = rep(1:nrow(perm), 2), x = rep(c(1, -1), each = nrow(perm)), dims = c(M + H, nrow(perm)))
	
			# updating E[Z], the allocation of cells on prototypes
			Z <- as.matrix(-1/2  * landscape[['beta']] * rdist(t(V), t(Theta))^2 - 1/2 * K * log(2 * pi) + 1/2 * K * log(landscape[['beta']]) + matrix(landscape[['pi.log']], M, H, byrow = TRUE))
			if (!is.null(W))
				Z[!W] <- -Inf
			Z <- t(softmax(t(Z)))# E[Z]
	
			BETA.blk <- kronecker(Diagonal(n = M), matrix(-landscape[['beta']], nrow = K, ncol = K))	# helper matrix for computing hessian
			iter <- 1
			while (iter <= max.iter){	# Laplace approximation
				V0 <- V
				PD.log <- V + B
				# compute the gradient of cells to update
				GR <- SN - exp(PD.log) - landscape[['beta']] * matrix(Matrix::rowSums(matrix(Matrix::colSums(cbind(V, Theta) %*% P), nrow = M, ncol = H) * Z), nrow = K, ncol = M, byrow = TRUE) - V
				# compute the inverse of covariance matrix using Sherman-Morrison formula
				PD.inv <- 1 / (exp(PD.log) + 1)
				w <- colSums(PD.inv)
				PD.INV <- Diagonal(x = -c(PD.inv))
				H.INV <- PD.INV - PD.INV %*% BETA.blk %*% PD.INV %*% kronecker(Diagonal(x = 1 / (1 + landscape[['beta']] * w)), Diagonal(n = K))
				# create a block diagonal matrix for the hession matrix and update V
				V <- V - matrix(H.INV %*% c(GR), K, M)	# E[V]
				V.var <- matrix(-diag(H.INV), K, M)	# diagonal of Var[v_m], for computing E[exp(V)]
				if (max.iter > 1 && norm(V - V0, '2') < 1e-7)
					break
				iter <- iter + 1
			}
			list(
				a = a, 
				V = V, 
				V.exp = exp(V + 1/2 * V.var),	# E[exp(V)]
				Z = Z
			)
		}else if (model == 'ccm'){
		
			B <- matrix(a, nrow = K, ncol = M, byrow = TRUE) + log(t(U) %*% matrix(1, N, M) + eps) # sum of a[m] + log(U[n, k] * W[n, m]) from n=1 to n=N
			iter <- 1
			while (iter <= max.iter){	# Laplace approximation
				V0 <- V
				PD.log <- V + B
				GR <- as.matrix(SN - exp(PD.log) - V)
				# compute the inverse using Sherman-Morrison formula
				PD.inv <- 1 / (exp(PD.log) + 1)
				H.INV <- Diagonal(x = -c(PD.inv))
				V <- V - matrix(H.INV %*% c(GR), K, M)	# E[V]
				V.var <- matrix(-diag(H.INV), K, M)	# diagonal of Var[v_m], for computing E[exp(V)]
				if (max.iter > 1 && norm(V - V0, '2') < 1e-7)
					break
				iter <- iter + 1
			}
			list(
				a = a, 
				V = V, 
				V.exp = exp(V + 1/2 * V.var)	# E[exp(V)]
			)
		}else if (model == 'ccr'){
			list(a = a, V.exp = exp(V))
		}
	}else if (model == 'gtm'){
		H <- landscape[['H.prototype']]	# number of total prototypes
		M <- ncol(V)	#	number of cells
		K <- nrow(V)
		# updating E[Z], the allocation of cells on prototypes
		Z <- as.matrix(-1/2  * landscape[['beta']] * rdist(t(V), t(landscape[['Theta.free']] %*% landscape[['S']]))^2 - 1/2 * K * log(2 * pi) + 1/2 * K * log(landscape[['beta']]) + matrix(landscape[['pi.log']], M, H, byrow = TRUE))
		if (!is.null(W))
			Z[!W] <- -Inf
		Z <- t(softmax(t(Z)))# E[Z]
		list(V = V, Z = Z)
	}
} # end of update.localvar


#' describe the prototype landscape
#' 
#' @param type type of the prototype landscape
#' @param n.circle prototypes per layer (S)
#' @param n.prototype the number of layers (R)
#' @param n.prev the number of layers of convolving prototypes (R - rho)
init.landscape <- function(type = 'temporal.convolving', K = 10, ...){

	param <- list(...)

	if (is.null(param[['alpha0']]))
		alpha0 <- 0.001
	else
		alpha0 <- param[['alpha0']]

	if (type == 'temporal.convolving'){

		time.points <- param[['time.points']]

		if (is.null(time.points) || time.points < 2)
			stop('at least two time points must be specified for landscape tcm')

		n.prototype <- param[['n.prototype']]
		n.circle <- param[['n.circle']]
		n.prev <- param[['n.prev']]
		
		if (is.null(param[['lambda']]))
			lambda <- 1
		else
			lambda <- param[['lambda']]

		H <- n.circle * n.prototype # number of prototypes per time point
		H.prototype <- H * time.points	# total number of prototypes

		# 2D coordinate of a generic prototype landscape
		Y0p <- as.matrix(expand.grid(angle = seq(0, 2 * pi, length.out = n.prototype + 1)[-1], r = seq(0, 1, length.out = n.circle + 1)[-1]))	# polar coord
		Y0 <- rbind(x = Y0p[, 'r'] * cos(Y0p[, 'angle']), y = Y0p[, 'r'] * sin(Y0p[, 'angle']))	# x-y coord

		L0.inv <- rdist(t(Y0))	# pairwise euclidean distance between prototypes
		L0.inv <- exp(-L0.inv^2 / (2 * lambda^2))	# squared exponential covariance function
		L0.inv <- as.matrix(nearPD(L0.inv)$mat)	# convert to a nearest PSD matrix
		L0 <- chol2inv(chol(L0.inv))	# inverse the corrected covariance matrix by Cholesky decomposition
		L.inv <- kronecker(diag(rep(1, time.points)), L0.inv)	# prototypes from every layer (time point)
		L.inv <- as.matrix(nearPD(L.inv)$mat)	# convert to a nearest PSD matrix
		L <- chol2inv(chol(L.inv))	# inverse the corrected covariance matrix by Cholesky decomposition

		TMC <- sparseMatrix(i = rep(1:time.points, each = H), j = 1:H.prototype, dims = c(time.points, H.prototype))	# time ~ prototypes
		MC <- NULL

		is.free0 <- Y0p[, 'r'] >  n.prev / n.circle	# the free prototypes at current layer
		Y0.nc <- Y0[, !is.free0]	# the 2D coordinate of the convolving prototypes for current layer
		Y0.nc <- Y0.nc * n.circle / n.prev	# scale to [-1, 1]
		S.conv <- softmax(-rdist(t(Y0), t(Y0.nc))^2)	# a stochastic matrix mapping from prototypes from layer t to non-scaffold prototypes from layer t + 1

		is.free <- c(rep(TRUE, H), rep(is.free0, time.points - 1))	# indicator for free prototypes
		H.free <- sum(is.free)
		S <- Diagonal(n = H.prototype)
    S <- S[is.free, ] # a H.free ~ H.prototype stochastic matrix mapping from free prototypes to ALL prototypes

		for (t in 2:time.points){
			hc <- TMC[t, ] 	# prototypes from current time point
			hp <- TMC[t - 1, ] 	# prototypes from previous time point
			hc.conv <- hc & !is.free # convolving prototypes from current time point
			S[, hc.conv] <- S[, hp] %*% S.conv
		}

		L <- S %*% L %*% t(S)
		
		# establish a global 2D coordinate of prototypes from all layers/time points
		Yp.prototype <- kronecker(rep(1, time.points), Y0p)
		colnames(Yp.prototype) <- c('angle', 'r')
		Yp.prototype[, 2] <- Yp.prototype[, 2] + rep(0:(time.points - 1), each = H)
		Y.prototype <- rbind(x = Yp.prototype[, 'r'] * cos(Yp.prototype[, 'angle']), y = Yp.prototype[, 'r'] * sin(Yp.prototype[, 'angle']))
		
		Theta.free <- matrix(rnorm(K * H.free), K, H.free)
		pi.log <- rep(digamma(alpha0) - digamma(alpha0 * H.prototype), H.prototype)	# E[log(pi)]

	}else if (type == 'plate'){

		if (is.null(param[['lambda']]))
			lambda <- 1
		else
			lambda <- param[['lambda']]

		n.prototype <- param[['n.prototype']]
		n.circle <- param[['n.circle']]
		H.free <- H.prototype <- H <- n.circle * n.prototype # number of prototypes 
		Yp.prototype <- as.matrix(expand.grid(angle = seq(0, 2 * pi, length.out = n.prototype + 1)[-1], r = seq(0, 1, length.out = n.circle + 1)[-1]))	# polar coord
		Y.prototype <- rbind(x = Yp.prototype[, 'r'] * cos(Yp.prototype[, 'angle']), y = Yp.prototype[, 'r'] * sin(Yp.prototype[, 'angle']))	# x-y coord
		L.inv <- rdist(t(Y.prototype))	# pairwise euclidean distance between prototypes
		L.inv <- exp(-L.inv^2 / (2 * lambda^2))	# squared exponential covariance function
		L.inv <- as.matrix(nearPD(L.inv)$mat)	# convert to a nearest PSD matrix
		L <- chol2inv(chol(L.inv))	# inverse the corrected covariance matrix by Cholesky decomposition
		S <- Diagonal(n = H.prototype)
		is.free <- rep(TRUE, H.prototype)
		n.prev <- NULL; TMC <- NULL
		MC <- sparseMatrix(i = 1:H.prototype, j = rep(1:n.circle, each = n.prototype), dims = c(H.prototype, n.circle))	# prototype ~ circle

		Theta.free <- matrix(rnorm(K * H.free), K, H.free)
		pi.log <- rep(digamma(alpha0) - digamma(alpha0 * H.prototype), H.prototype)	# E[log(pi)]

	}else
		stop(sprintf('unknown landscape type: %s', type))

	structure(list(
		type = type,
		K = K,
		H = H, 
		H.prototype = H.prototype, 
		n.prototype = n.prototype, 
		n.circle = n.circle, 
		n.prev = n.prev, 
		lambda = lambda, 
		TMC = TMC, 
		MC = MC,
		L = L, 
		S = S, 
		H.free = H.free, 
		Theta.free = Theta.free, 
		Y.prototype = Y.prototype,
		Yp.prototype = Yp.prototype,
		alpha0 = alpha0,
		is.free = is.free,
		pi.log = pi.log,
		beta = 1
	), class = 'landscape')

} # end of init.landscape 


#' Correlated Cell Model
#'
ccm <- function(X, U = NULL, V = NULL, a = NULL, optimization.method = 'batch', c0 = 0.001, d0 = 0.001, max.iter = 100, verbose = TRUE, batch.size = 100, decay.rate = 1, forgetting.rate = 0.75, mc.cores = 1, max.size.per.batch = 1000){

	eps <- .Machine$double.eps

	# the number of iterations for Newton step for optimizing V and V.var
	max.iter2 <- switch(optimization.method, 'batch' = 50, 'stochastic' = 1)

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- nrow(V)	# number of metagenes

	# splitting the samples into batches
	if (M <= mc.cores * max.size.per.batch){	
		n.batch <- mc.cores
	}else{
		n.batch <- ceiling(M / max.size.per.batch)
	}
	groups <- sample(1:n.batch, M, replace = TRUE)
	size <- table(factor(groups, 1:n.batch))	# number of samples per core

	# initialize the local variables
	localvar <- mclapply(1:n.batch, function(b){
		m <- groups == b
		V.var <- matrix(1, nrow = K, ncol = size[b])	# Var[V]
		V.exp <- exp(V[, m] + 1/2 * V.var)	# E[exp(V)]
		list(a = a[m], V = V[, m], V.exp = V.exp)
	}, mc.cores = mc.cores)

	iter <- 1
	optval <- NULL
	while (iter <= max.iter){

		Up <- U

		# updating all local variables: E[V], E[exp(V)] and a
		localvar <- mclapply(1:n.batch, function(b){
			m <- groups == b
			update.localvar(
				X = X[, m, drop = FALSE],
				U = U,
				V = localvar[[b]]$V, 
				V.exp = localvar[[b]]$V.exp,
				model = 'ccm', max.iter = max.iter2
			)
		}, mc.cores = mc.cores)

		rho <- switch(optimization.method,
			'stochastic' = learning.rate(iter, decay.rate, forgetting.rate),
			'batch' = 1
		)

		if (optimization.method == 'stochastic'){
			# SVI: sampling a subset of samples(cells) for updating global variables
			m <- lapply(1:n.batch, function(b) sample(1:size[b], min(max(1, round(batch.size / n.batch)), size[b])))
		}else{
			# batch optimization: use all the samples
			m <- lapply(1:n.batch, function(b) 1:size[b])
		}

		U <- (1 - rho) * U + rho * Reduce('+', mclapply(1:n.batch, function(b){
			c <- U * (X[, which(groups == b)[m[[b]]]] / (U %*% localvar[[b]]$V.exp[, m[[b]], drop = FALSE] + eps)) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE]) + c0
			d <- matrix(1, N, length(m[[b]])) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE] %*% diag(exp(localvar[[b]]$a[m[[b]]]))) + d0
			U.log <- digamma(c) - log(d + eps)	# E[log(U)]
			exp(U.log)
		}, mc.cores = mc.cores)) / n.batch

		if (iter == 1 || iter %% 10 == 0){
			J <- norm(U - Up, '2')
			cat(sprintf('[%s] ccm | iter=%5.d | M=%5.d | n.batch=%3.d | dJ=%7.3e\n', Sys.time(), iter, M, n.batch, J))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), J = J))
			if (J < 0.01)
				break
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	a <- rep(0, M)
	V <- matrix(0, K, M)
	for (b in 1:mc.cores){
		V[, groups == b] <- localvar[[b]]$V
		a[groups == b] <- localvar[[b]]$a
	}
	list(N = N, M = M, K = K, U = U, V = V, a = a, optval = optval)

} # end of ccm


#' The generative topographic model
#'
gtm <- function(V, landscape = NULL, CT = NULL, optimization.method = 'batch', min.beta = 0, max.iter = 100, verbose = TRUE, batch.size = 100, decay.rate = 1, forgetting.rate = 0.75, mc.cores = 1, max.size.per.batch = 1000){

	eps <- .Machine$double.eps
	M <- ncol(V)	# number of cells
	K <- nrow(V)	# number of metagenes

	# the number of iterations for Newton step for optimizing V and V.var
	max.iter2 <- switch(optimization.method, 'batch' = 50, 'stochastic' = 1)
	
	if (!is.null(landscape)){
		Theta <- landscape[['Theta.free']] %*% landscape[['S']]	# all prototypes
		if (landscape[['type']] == 'temporal.convolving'){
			W <- CT %*% landscape[['TMC']]	# a binary cell ~ prototype assignment matrix
		}else if (landscape[['type']] %in% c('som', 'plate')){
			W <- matrix(TRUE, M, landscape[['H.prototype']])
		}
	}else
		stop('landscape must be specified')


	# splitting the samples into batches
	if (M <= mc.cores * max.size.per.batch){	
		n.batch <- mc.cores
	}else{
		n.batch <- ceiling(M / max.size.per.batch)
	}
	groups <- sample(1:n.batch, M, replace = TRUE)
	size <- table(factor(groups, 1:n.batch))	# number of samples per core
	landscape[['pi.log']] <- rep(digamma(landscape[['alpha0']]) - digamma(landscape[['alpha0']]* landscape[['H.prototype']]), landscape[['H.prototype']])

	iter <- 1
	optval <- NULL
	while (iter <= max.iter){

		Theta.free.p <- landscape[['Theta.free']]

		# updating all local variables: E[Z]
		localvar <- mclapply(1:n.batch, function(b){
			m <- groups == b
			update.localvar(
				V = V[, m, drop = FALSE],
				W = W[m, , drop = FALSE],
				landscape = landscape,
				model = 'gtm'
			)
		}, mc.cores = mc.cores)

		rho <- switch(optimization.method,
			'stochastic' = learning.rate(iter, decay.rate, forgetting.rate),
			'batch' = 1
		)

		if (optimization.method == 'stochastic'){
			# SVI: sampling a subset of samples(cells) for updating global variables
			m <- lapply(1:n.batch, function(b) sample(1:size[b], min(max(1, round(batch.size / n.batch)), size[b])))
		}else{
			# batch optimization: use all the samples
			m <- lapply(1:n.batch, function(b) 1:size[b])
		}

		# updating free prototypes Theta
		SZS <- landscape[['S']] %*% Diagonal(x = landscape[['beta']] * Reduce('+', mclapply(1:n.batch, function(b) colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]), mc.cores = mc.cores))) %*% t(landscape[['S']])
		SZSL.inv <- chol2inv(chol(as.matrix(SZS + landscape[['L']])))
		landscape[['Theta.free']] <- (1 - rho) * landscape[['Theta.free']] + rho * landscape[['beta']] * Reduce('+', mclapply(1:n.batch, function(b){
			localvar[[b]]$V[, m[[b]], drop = FALSE] %*% localvar[[b]]$Z[m[[b]], , drop = FALSE] %*% t(landscape[['S']]) %*% SZSL.inv
		}, mc.cores = mc.cores))
		Theta <- landscape[['Theta.free']] %*% landscape[['S']]	# all prototypes

		# updating beta
		landscape[['beta']] <- (1 - rho) * landscape[['beta']] + rho * Reduce('+', mclapply(1:n.batch, function(b){
			(K * sum(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + 0.001) / (sum(localvar[[b]]$Z[m[[b]], , drop = FALSE] * rdist(t(localvar[[b]]$V[, m[[b]], drop = FALSE]), t(Theta))^2) + 0.001)
		}, mc.cores = mc.cores)) / n.batch
		landscape[['beta']] <- max(landscape[['beta']], min.beta)

		# updating pi.log
		kappa <- Reduce('+', mclapply(1:n.batch, function(b){
			colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + landscape[['alpha0']]
		}, mc.cores = mc.cores)) 
		landscape[['pi.log']] <- (1 - rho) * landscape[['pi.log']] + rho * (digamma(kappa) - digamma(sum(kappa))) # E[log(pi)]
		
		if (iter == 1 || iter %% 10 == 0){
			J <- norm(landscape[['Theta.free']] - Theta.free.p, '2')
			n.cluster <- length(unique(unlist(mclapply(1:n.batch, function(b) max.col(localvar[[b]]$Z), mc.cores = mc.cores))))
			cat(sprintf('[%s] gtm | iter=%5.d | M=%5.d | n.batch=%3.d | dJ=%7.3e | # non-empty prototypes=%4.d | beta=%7.3f\n', Sys.time(), iter, M, n.batch, J, n.cluster, landscape[['beta']]))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), n.cluster = n.cluster, beta = landscape[['beta']], J = J))
			if (J < 0.01)
				break
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	Z <- matrix(0, M, landscape[['H.prototype']])
	for (b in 1:mc.cores){ Z[groups == b, ] <- localvar[[b]]$Z
	}

	structure(list(M = M, K = K, Z = Z, CT = CT, optval = optval, landscape = landscape), class = 'gtm')

} # end of gtm


#' log density of Poission distribution
#'
ldpois <- function(k, mu){

		mu <- mu + .Machine$double.eps

	if (any(k <= -1))
				stop('k must be greater than -1')
		
		k * log(mu) - mu - lgamma(k + 1)

} # end of ldpois

plot.gtm <- function(x, ...) plot.tcm(x, ...)

print.landscape <- function(x, ...){

	cat(sprintf('prototype landscape: %s\n', x[['type']]))
	cat(sprintf('number of metagenes(K): %d\n', x[['K']]))
	cat(sprintf('number of total prototypes: %d\n', x[['H.prototype']]))
	cat(sprintf('number of free prototypes: %d\n', x[['H.free']]))
	cat(sprintf('lambda: %.3e\n', x[['lambda']]))

	if (x[['type']] == 'temporal.convolving'){
		cat(sprintf('number of time points: %d\n', ncol(x[['TMC']])))
		cat(sprintf('number of circle per time point: %d\n', x[['n.circle']]))
		cat(sprintf('number of prototypes per circle: %d\n', x[['n.prototype']]))
		cat(sprintf('number of prototypes for each time point: %d\n', x[['H']]))
		cat(sprintf('number of circle(s) mapped from the previous time point: %d\n', x[['n.prev']]))
	}else if (x[['type']] == 'plate'){
		cat(sprintf('number of circles: %d\n', x[['n.circle']]))
		cat(sprintf('number of prototypes per circle: %d\n', x[['n.prototype']]))
	}
} # end of print.landscape


#' Correlated Cell Regression
#'
ccr <- function(X, U = NULL, V = NULL, a = NULL, optimization.method = 'batch', c0 = 0.001, d0 = 0.001, max.iter = 100, verbose = TRUE, batch.size = 100, decay.rate = 1, forgetting.rate = 0.75, mc.cores = 1, max.size.per.batch = 1000){

	eps <- .Machine$double.eps

	# the number of iterations for Newton step for optimizing V and V.var
	max.iter2 <- switch(optimization.method, 'batch' = 50, 'stochastic' = 1)

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- nrow(V)	# number of metagenes

	# splitting the samples into batches
	if (M <= mc.cores * max.size.per.batch){	
		n.batch <- mc.cores
	}else{
		n.batch <- ceiling(M / max.size.per.batch)
	}
	groups <- sample(1:n.batch, M, replace = TRUE)
	size <- table(factor(groups, 1:n.batch))	# number of samples per core

	# initialize the local variables
	localvar <- mclapply(1:n.batch, function(b) list(a = a[groups == b]), mc.cores = mc.cores)
	if (is.null(U))
		U <- matrix(runif(N * K), N, K)

	iter <- 1
	optval <- NULL
	Jp <- 0
	while (iter <= max.iter){

		# updating all local variables: E[V], E[exp(V)] and a
		localvar <- mclapply(1:n.batch, function(b){
			m <- groups == b
			update.localvar(
				X = X[, m, drop = FALSE],
				U = U,
				V = V[, m, drop = FALSE],
				model = 'ccr', max.iter = max.iter2
			)
		}, mc.cores = mc.cores)

		rho <- switch(optimization.method,
			'stochastic' = learning.rate(iter, decay.rate, forgetting.rate),
			'batch' = 1
		)

		if (optimization.method == 'stochastic'){
			# SVI: sampling a subset of samples(cells) for updating global variables
			m <- lapply(1:n.batch, function(b) sample(1:size[b], min(max(1, round(batch.size / n.batch)), size[b])))
		}else{
			# batch optimization: use all the samples
			m <- lapply(1:n.batch, function(b) 1:size[b])
		}

		U <- (1 - rho) * U + rho * Reduce('+', mclapply(1:n.batch, function(b){
			c <- U * (X[, which(groups == b)[m[[b]]]] / (U %*% localvar[[b]]$V.exp[, m[[b]], drop = FALSE] + eps)) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE]) + c0
			d <- matrix(1, N, length(m[[b]])) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE] %*% diag(exp(localvar[[b]]$a[m[[b]]]))) + d0
			U.log <- digamma(c) - log(d + eps)	# E[log(U)]
			exp(U.log)
		}, mc.cores = mc.cores)) / n.batch

		if (iter == 1 || iter %% 10 == 0){
			J <- sum(unlist(mclapply(1:n.batch, function(b){
				sum(ldpois(X[, groups == b], U %*% localvar[[b]]$V.exp %*% Diagonal(x = exp(localvar[[b]]$a))))
			}, mc.cores = mc.cores)))
			cat(sprintf('[%s] ccr | iter=%5.d | M=%5.d | n.batch=%3.d | dJ=%7.3e\n', Sys.time(), iter, M, n.batch, abs(J - Jp)))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), J = J))
			if (abs(J - Jp) < 0.01)
				break
			Jp <- J
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	a <- rep(0, M)
	for (b in 1:mc.cores){
		a[groups == b] <- localvar[[b]]$a
	}
	list(N = N, M = M, K = K, U = U, a = a, optval = optval)

} # end of ccr

