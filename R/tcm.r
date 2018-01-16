library(irlba)	# irlba
library(fields)	# rdist
library(gtools)	# ddirichlet
library(FNN)	# knnx.index
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
tcm <- function(X, K = 10, time.table = NULL, landscape = NULL, optimization.method = 'stochastic', batch.size = 2000, max.iter = 100, mc.cores = 1){

	N <- nrow(X)
	M <- ncol(X)
	eps <- .Machine$double.eps

	if (is.null(time.table) && landscape[['type']] == 'temporal.convolving')
		stop('time.table cannot be NULL if landscape$type is temporal.convolving')
	
	if (is.null(landscape))
		stop('landscape cannot be NULL')

	CT <- as(as(time.table, 'matrix'), 'ngCMatrix') # cell ~ time
	time <- colnames(time.table)

	cat('-------------------------------------------------------------------------------------------------\n')
	cat('topographic cell map\n')
	cat('-------------------------------------------------------------------------------------------------\n')
	cat(sprintf('[%s] number of input genes(N): %d\n', Sys.time(), N))
	cat(sprintf('[%s] number of input cells(M): %d\n', Sys.time(), M))
	cat(sprintf('[%s] number of metagenes(K): %d\n', Sys.time(), K))
	cat(sprintf('[%s] optimization method: %s\n', Sys.time(), optimization.method))
	cat(sprintf('[%s] number of cores: %d\n', Sys.time(), mc.cores))
	if (!is.null(time.table))
		cat(sprintf('[%s] number of time points: %d\n', Sys.time(), ncol(CT)))

	cat(sprintf('[%s] initializing metagene coefficients and basis\n', Sys.time()))
	V <- scale(t(irlba(scale(log(X + 1)), nu = 1, nv = K)$v))
	a <- rep(0, M)	# the cell effect
	U <- matrix(runif(N * K), N, K)

	if (landscape[['type']] == 'temporal.convolving'){	# initialization for landscape type = temporal.convolving
		for (t in ncol(CT):2){
			m <- Matrix::rowSums(CT[, t:ncol(CT), drop = FALSE]) > 0	#  cells for training TCM
			cat(sprintf('[%s] initializing %d cells from time point(s) %d-%d\n', Sys.time(), sum(m), t, ncol(CT)))
			mf <- tcm.core(X = X[, m], U = U, V = V[, m], a = a[m], landscape = landscape, CT = CT[m, ], optimization.method = optimization.method, update.beta = TRUE, batch.size = batch.size, max.iter = 50, mc.cores = mc.cores)
			U <- mf$U
			V[, m] <- mf$V
			a[m] <- mf$a
			landscape$Theta.free <- mf$Theta.free
		}
		cat(sprintf('[%s] fitting TCM\n', Sys.time()))
		mf <- tcm.core(X = X, U = U, V = V, a = a, landscape = landscape, CT = CT, optimization.method = optimization.method, batch.size = batch.size, max.iter = max.iter, mc.cores = mc.cores)
		mf$CT <- CT
		mf$time <- time
	}

	mf$landscape <- landscape
	class(mf) <- 'tcm'
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
#' @param test.size number of random samples for evaluating objective function (default: 100)
#' @param mc.cores number of CPU cores (default: 1)
#' @param max.size.per.batch maximum samples per batch
#'
tcm.core <- function(X, U = NULL, V = NULL, a = NULL, landscape = NULL, CT = NULL, C = NULL, optimization.method = 'batch', beta = 1, update.beta = TRUE, c0 = 0.001, d0 = 0.001, max.iter = 100, verbose = TRUE, batch.size = 100, decay.rate = 1, forgetting.rate = 0.75, test.size = 100, mc.cores = 1, max.size.per.batch = 1000){

	eps <- .Machine$double.eps

	# the number of iterations for Newton step for optimizing V and V.var
	max.iter2 <- switch(optimization.method, 'batch' = 50, 'stochastic' = 1)
	
	if (!is.null(landscape)){
		Theta.free <- landscape$Theta.free
		L <- landscape$L
		S <- landscape$S
		Theta <- Theta.free %*% S	# all prototypes
		W <- CT %*% landscape$TMC	# a binary cell ~ prototype assignment matrix
	}else
		stop('landscape must be specified')

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- ncol(U)	# number of metagenes

	if (is.null(C))
		C <- sparseMatrix(i = 1:M, j = rep(1, M), dims = c(M, 1))

	if (is.null(a))
		a <- rep(0, M)	# the cell effect
	
	nc <- ncol(C)	# number of predefined cell groups
	H <- ncol(S)	# number of prototypes

	Alpha0 <- matrix(landscape$alpha0, nc, H)
	Pi.log <- digamma(Alpha0) - matrix(digamma(Matrix::rowSums(Alpha0)), nc, H)	# initialize E[log(pi)], nc by H

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

		Theta.free.p <- Theta.free

		# updating all local variables: E[V], E[exp(V)], a and E[Z]
		localvar <- mclapply(1:n.batch, function(b){
			m <- groups == b
			update.localvar(
				X = X[, m, drop = FALSE],
				U = U,
				V = localvar[[b]]$V, 
				V.exp = localvar[[b]]$V.exp,
				Theta = Theta,
				beta = beta,
				C = C[m, ],
				Pi.log = Pi.log,
				W = W[m, ],
				max.iter = max.iter2
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
		SZS <- S %*% Diagonal(x = beta * Reduce('+', mclapply(1:n.batch, function(b) colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]), mc.cores = mc.cores))) %*% t(S)
		tryCatch({
			SZSL.inv <- chol2inv(chol(as.matrix(SZS + L)))
		}, error = function(e) browser())	
		Theta.free <- (1 - rho) * Theta.free + rho * beta * Reduce('+', mclapply(1:n.batch, function(b){
			localvar[[b]]$V[, m[[b]], drop = FALSE] %*% localvar[[b]]$Z[m[[b]], , drop = FALSE] %*% t(S) %*% SZSL.inv
		}, mc.cores = mc.cores))
		Theta <- Theta.free %*% S	# mapping free prototypes to all prototypes

		U <- (1 - rho) * U + rho * Reduce('+', mclapply(1:n.batch, function(b){
			c <- U * (X[, which(groups == b)[m[[b]]]] / (U %*% localvar[[b]]$V.exp[, m[[b]], drop = FALSE] + eps)) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE]) + c0
			d <- matrix(1, N, length(m[[b]])) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE] %*% diag(exp(localvar[[b]]$a[m[[b]]]))) + d0
			U.log <- digamma(c) - log(d + eps)	# E[log(U)]
			exp(U.log)
		}, mc.cores = mc.cores)) / n.batch

		if (update.beta){
			# updating beta
			beta <- (1 - rho) * beta + rho * Reduce('+', mclapply(1:n.batch, function(b){
				(K * sum(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + 0.001) / (sum(localvar[[b]]$Z[m[[b]], , drop = FALSE] * rdist(t(localvar[[b]]$V[, m[[b]], drop = FALSE]), t(Theta))^2) + 0.001)
			}, mc.cores = mc.cores)) / n.batch
		}

		# updating Kappa/Pi.log
		Kappa <- Reduce('+', mclapply(1:n.batch, function(b){
			as.matrix(t(C[m[[b]], , drop = FALSE]) %*% localvar[[b]]$Z[m[[b]], , drop = FALSE] + Alpha0)
		}, mc.cores = mc.cores)) 
		Pi.log <- (1 - rho) * Pi.log + rho * (digamma(Kappa) - matrix(digamma(rowSums(Kappa)), nc, H)) # E[log(pi)]

		if (iter == 1 || iter %% 10 == 0){
			n.cluster <- length(unique(unlist(mclapply(1:n.batch, function(b) max.col(localvar[[b]]$Z), mc.cores = mc.cores))))
			cat(sprintf('[%s] iter=%5.d | M=%d | n.batch=%d | # non-empty prototypes=%4.d | beta=%7.3f | d||Theta.free||=%7.3f\n', Sys.time(), iter, M, n.batch, n.cluster, beta, norm(Theta.free - Theta.free.p, '2')))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), n.cluster = n.cluster, beta = beta))
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
	list(H = H, N = N, M = M, K = K, U = U, V = V, a = a, Z = Z, Theta.free = Theta.free, C = C, Alpha0 = Alpha0, optval = optval)

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

	if (mf$landscape$type == 'temporal.convolving'){

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

		for (t in 1:ncol(CT))
			draw.circle(0, 0, t, nv = 100, border = 'black', lty = 2, lwd = 2)

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
update.localvar <- function(X, U, V = NULL, V.exp = NULL, Theta = NULL, beta = NULL, C = NULL, Pi.log = NULL, W = NULL, max.iter = 1){

	eps <- .Machine$double.eps
	K <- ncol(U)
	N <- nrow(U)
	M <- ncol(X)

	if (is.null(V))
		V <- matrix(rnorm(K * M), K, M)

	if (is.null(V.exp))
		V.exp <- exp(V)

	W <- as.matrix(W)

	SN <- as.matrix(V.exp * ( t(U) %*% ( X / (U %*% V.exp + eps)) ))
	a <- log(colSums(SN)) - log(colSums(U %*% V.exp) + eps)	# updating the cell effect (a)
	B <- matrix(a, nrow = K, ncol = M, byrow = TRUE) + log(t(U) %*% matrix(1, N, M) + eps) # sum of a[m] + log(U[n, k] * W[n, m]) from n=1 to n=N

	H <- ncol(Theta)	# number of prototypes
	# construct a permutation matrix to compute \sum_{m=1}^M \sum_{h=1}^H ||\Theta_h - V_m||^2
	perm <- expand.grid(1:M, (M + 1):(M + H))
	P <- sparseMatrix(i = c(perm[, 1], perm[, 2]), j = rep(1:nrow(perm), 2), x = rep(c(1, -1), each = nrow(perm)), dims = c(M + H, nrow(perm)))

	# updating E[Z], the allocation of cells on prototypes
	Z <- as.matrix(-1/2  * beta * rdist(t(V), t(Theta))^2 - 1/2 * K * log(2 * pi) + 1/2 * K * log(beta) + C %*% Pi.log)
	if (!is.null(W))
		Z[!W] <- -Inf
	Z <- t(softmax(t(Z)))# E[Z]

	BETA.blk <- kronecker(Diagonal(n = M), matrix(-beta, nrow = K, ncol = K))	# helper matrix for computing hessian
	iter <- 1
	while (iter <= max.iter){	# Laplace approximation
		V0 <- V
		PD.log <- V + B
		# compute the gradient of cells to update
		GR <- SN - exp(PD.log) - beta * matrix(Matrix::rowSums(matrix(Matrix::colSums(cbind(V, Theta) %*% P), nrow = M, ncol = H) * Z), nrow = K, ncol = M, byrow = TRUE) - V
		# compute the inverse of covariance matrix using Sherman-Morrison formula
		PD.inv <- 1 / (exp(PD.log) + 1)
		w <- colSums(PD.inv)
		PD.INV <- Diagonal(x = -c(PD.inv))
		H.INV <- PD.INV - PD.INV %*% BETA.blk %*% PD.INV %*% kronecker(Diagonal(x = 1 / (1 + beta * w)), Diagonal(n = K))
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
} # end of update.localvar


#' describe the prototype landscape
#' 
#' @param type type of the prototype landscape
#' @param n.circle prototypes per layer (S)
#' @param n.prototype the number of layers (R)
#' @param n.prev the number of layers of convolving prototypes (R - rho)
init.landscape <- function(type = 'temporal.convolving', ...){

	param <- list(...)
	
	if (type == 'temporal.convolving'){

		time.points <- param[['time.points']]

		if (is.null(time.points) || time.points < 2)
			stop('at least two time points must be specified for landscape tcm')

		K <- param[['K']]
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
		alpha0 <- 0.001

		cat(sprintf('prototype landscape: %s\n', type))
		cat(sprintf('number of metagenes(K): %d\n', K))
		cat(sprintf('number of time points: %d\n', time.points))
		cat(sprintf('number of circle per time point: %d\n', n.circle))
		cat(sprintf('number of prototypes per circle: %d\n', n.prototype))
		cat(sprintf('number of prototypes for each time point: %d\n', H))
		cat(sprintf('number of circle(s) mapped from the previous time point: %d\n', n.prev))
		cat(sprintf('number of total prototypes: %d\n', H.prototype))

	}else
		stop(sprintf('unknown landscape type: %s', type))

	list(
		type = type,
		H = H, 
		H.prototype = H.prototype, 
		n.prototype = n.prototype, 
		n.circle = n.circle, 
		n.prev = n.prev, 
		lambda = lambda, 
		TMC = TMC, 
		L = L, 
		S = S, 
		H.free = H.free, 
		Theta.free = Theta.free, 
		Y.prototype = Y.prototype,
		Yp.prototype = Yp.prototype,
		alpha0 = alpha0,
		is.free = is.free
	)

} # end of landscape 

