#' TCM: A package for visualizing temporal scRNA-seq data
#'
#' @import Matrix
#' @import SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom irlba irlba
#' @importFrom fields rdist
#' @importFrom gtools ddirichlet rdirichlet
#' @importFrom FNN knnx.index
#' @importFrom cluster pam
#' @importFrom MASS mvrnorm
#' @importFrom igraph get.shortest.paths graph.adjacency
#' @importFrom S4Vectors DataFrame
#' @importFrom gplots colorpanel
#' @importFrom plotrix draw.circle
#' @importFrom parallel mclapply
#' @docType package
#' @name tcm
NULL



#' add.paths
#'
#' A generic function for adding developmental trajectory onto the 
#' the prototype landscape
#'
#' @param x a tcm object
#' @param ... additional arguments
#'
#' @export
#'
#' @author Wuming Gong, \email{gongx030@umn.edu}
#'
add.paths <- function(x, ...) UseMethod('add.paths', x)


#' as.igraph
#'
#' A generic function for converting prototype landscape to a igraph object
#'
#' @param x an object
#' @param ... additional arguments
#'
#' @author Wuming Gong, \email{gongx030@umn.edu}
#'
as.igraph <- function(x, ...) UseMethod('as.igraph', x)


#' sim.rnaseq.ts
#'
#' Simulating time series single cell RNA-seq data
#' 
#' @param N number of genes
#' @param M number of cells
#' @param ls the prototype landscape (default: a "plate" with 100 circles and the number of prototypes per circle of 10)
#' @param lambda the argument of a exponential decay model for introducing the dropout noise (default: 0.25)
#' @param alpha0 the argument for a Dirichlet distribution for sampling metagene basis (default: 0.1)
#' @param library.size total number of reads per cell (default: c(1e3, 1e5))
#' @param n.lineage number of simulated lineages (default: 5)
#' @param type the type of differentiation models (default: 'sequential')
#' @param n.time.points the number of simulated time points (default: 5)
#' @param ... additional arguments
#' 
#' @return A SummarizedExperiment object of simulated temporal scRNA-seq data
#'
#' @export
#'
#' @seealso \code{\link{landscape}}
#'
#' @author Wuming Gong, \email{gongx030@umn.edu}
#' 
#' @examples
#' # simulate a simple temporal scRNA-seq data with 2,000 genes, 500 cells and five different lineages.  
#' # The single cell data are sampled across five time points following a sequentail differentiation model. 
#' set.seed(122)
#' sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, n.time.points = 5)
#'
sim.rnaseq.ts <- function(N = 2000, M = 500, ls, lambda.dropout = 0.25, alpha0 = 0.1, library.size = c(1e3, 1e5), n.lineage = 5, type = 'sequential', n.time.points = 5, ...){

	param <- list(...)

	if (missing(ls)){
		ls <- landscape(type = 'plate', K = 15, n.prototype = 10, n.circle = 100)
	}

	# selecting the leave prototypes
	border.prototypes <- which(ls$MC[, ls$n.circle])	# the prototypes at the outer region
	leaves <- round(seq(min(border.prototypes), max(border.prototypes), length.out = n.lineage + 1))[-1]

	# find the shortest paths from the origin to the leaves
	ls$paths <- lapply(leaves, function(leaf) as.vector(get.shortest.paths(ls$graph, from = 1, to = leaf)$vpath[[1]]))

	# a matrix for the active prototypes for each leave
	LN <- sparseMatrix(i =  unlist(ls$paths), j = rep(1:n.lineage, sapply(ls$paths, length)), dims = c(ls$H.prototype, n.lineage))
	ls$is.active <- Matrix::rowSums(LN) > 0

	if (type == 'sequential'){
		starts <- round(seq(1, ls$n.circle, length.out = n.time.points + 1))[-(n.time.points + 1)]	# start of time interval
		ends <- c(starts[-1] - 1, ls$n.circle)
		P <- do.call('cbind', lapply(1:n.time.points, function(i) 1:ls$n.circle >= starts[i] & 1:ls$n.circle <= ends[i]))	# n.circle ~ n.time.points
		P <- as.matrix(P %*% Diagonal(x = 1 / Matrix::colSums(P))) # column stochastic matrix of the distribution of time points on each circle
	}else if (type == 'delayed'){
		starts <- rep(1, n.time.points)
		ends <- round(seq(1, ls$n.circle, length.out = n.time.points + 1))[-1]
		P <- do.call('cbind', lapply(1:n.time.points, function(i) 1:ls$n.circle >= starts[i] & 1:ls$n.circle <= ends[i]))	# n.circle ~ n.time.points
		P <- as.matrix(P %*% Diagonal(x = 1 / Matrix::colSums(P))) # column stochastic matrix of the distribution of time points on each circle
	}else if (type == 'forward'){
		starts <- round(seq(1, ls$n.circle, length.out = n.time.points + 1))[-(n.time.points + 1)]	# start of time interval
		ends <- rep(ls$n.circle, n.time.points)
		P <- do.call('cbind', lapply(1:n.time.points, function(i) 1:ls$n.circle >= starts[i] & 1:ls$n.circle <= ends[i]))	# n.circle ~ n.time.points
		P <- as.matrix(P %*% Diagonal(x = 1 / Matrix::colSums(P))) # column stochastic matrix of the distribution of time points on each circle
	}else if (type == 'uniform'){
		P <- matrix(1 / ls$n.circle, ls$n.circle, n.time.points)
	}else if (type == 'mixed'){
		P <- do.call('cbind', lapply(1:n.time.points, function(i){
			a <- runif(1, min = 0, max = 5)
			b <- runif(1, min = 0, max = 5)
			breaks <- seq(0, 1, length.out = ls$n.circle + 1)
			breaks[1] <- breaks[1] - 1
			breaks[length(breaks)] <- breaks[length(breaks)] + 1
			table(cut(rbeta(10000, a, b), breaks = breaks)) / 10000
		}))
	}else
		stop(sprintf('unknown type: %s', type))

	GBP <- as.matrix(ls$MC %*% P) / ls$n.prototype	# metacell ~ time, column stochastic matrix of sampling probability

	z <- rep(NA, M)				# metacell assignment for each cell
	lineage <- rep(NA, M)	# lineage assignment for each cell
	time <- rep(NA, M)		# time assignment for each cell
	for (m in 1:M){
		lineage[m] <- sample(1:n.lineage, 1)	# randomly assign a lineage
		time[m] <- sample(1:n.time.points, 1)	# randomly assign a time point
		# randomly sample a prototype along current lineage path, according to the probability of time labels on the path
		z[m] <- sample(which(LN[, lineage[m]]), 1, prob = GBP[LN[, lineage[m]], time[m]])
	}
	V <- as.matrix(ls$Theta.free %*% ls$S[, z])	# metacell coefficient for each cell
	V.exp <- exp(V) %*% diag(1 / colSums(exp(V)))
	libsize <- runif(M, library.size[1], library.size[2])	# sampling the library size for each cell
	U <- t(rdirichlet(ls$K, alpha = rep(alpha0, N)))	# sampling the metagene basis
	Mu <- do.call('cbind', lapply(1:M, function(m) rmultinom(1, libsize[m], prob = U %*% V.exp[, m, drop = FALSE])))	# sampling the number of reads
	prob <- exp(-lambda.dropout * log(Mu + 1))	# introduce the dropout noise by a exponential decay model
	D <- matrix(rbinom(N * M, 1, prob), nrow = N, ncol = M) == 1
	X <- Mu
	X[D] <- 0
	rownames(X) <- rownames(Mu) <- sprintf('G%d', 1:N)

	SummarizedExperiment(
		assays <- list(count = X, Mu = as.matrix(Mu), D = D),
		rowData = DataFrame(U = I(U)),
		colData = DataFrame(V = I(t(V)), z = z, lineage = lineage, circle = max.col(ls$MC)[z], time = I(time), time.table = I(table(1:M, factor(time)))),
		metadata = list(N = N, M = M, landscape = ls, n.lineage = n.lineage, lambda.dropout = lambda.dropout, type = type, n.time.points = n.time.points, P = P)
	)

} # end of sim.rnaseq.ts


#' tcm
#' 
#' Visualizing temporal scRNA-seq by topographic cell map 
#' 
#' @param X a read count matrix where each row represents a gene and each column represents a cell
#' @param time.table a cell by time point table indicating the source of each cell
#' @param ls the prototype landscape for TCM 
#'				(default: landscape(type = 'temporal.convolving', time.points = ncol(CT), K = 15, n.prototype = 15, n.circle = 10, n.prev = 3))
#' @param init: initialization parameters
#' @param	control: control parameters
#' 
#' @return a tcm object
#'
#' @export
#
#' @author Wuming Gong, \email{gongx030@umn.edu}
#
#' @examples
#'
#' library(tcm)
#'
#' # simulate a simple temporal scRNA-seq data with 2,000 genes, 500 cells and five different lineages.  
#' # The single cell data are sampled across five time points following a sequentail differentiation model. 
#' set.seed(122)
#' sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, n.time.points = 5)
#' X <- assays(sim)$count
#' time.table <- colData(sim)$time.table
#' mf <- tcm(X, time.table = time.table)
#' bg.cell <- rainbow(5)[colData(sim)$lineage]
#' dev.new(height = 10, width = 12)
#' par(mar = c(5, 5, 5, 15))
#' plot(mf, pch = 21, bg = bg.cell, cex = 1.5)
#' 
#' legend(par('usr')[2], par('usr')[4], 1:sim$n.lineage, bty = 'n', xpd = NA, pt.bg = bg.lineage, pch = 21, col = 'black', cex = 1.75)
#' 
tcm <- function(X, time.table, ls, init = NULL, control = NULL){

	N <- nrow(X)
	M <- ncol(X)
	eps <- .Machine$double.eps
	control <- set.control(control)

	if (missing(time.table))
		CT <- as(matrix(TRUE, M, 1), 'ngCMatrix')
	else	
		CT <- as(as(time.table, 'matrix'), 'ngCMatrix') # cell ~ time

	if (M != nrow(CT))
		stop('ncol(X) must be equal to nrow(time.table)')

	if (missing(ls)){
		cat(sprintf('[%s] landscape(ls) is missing\n', Sys.time()))
		cat(sprintf('[%s] using a temporal convolving landscape with K=15, n.prototype=15, n.circle=10, n.prev=3\n', Sys.time()))
		ls <- landscape(type = 'temporal.convolving', time.points = ncol(CT), K = 15, n.prototype = 15, n.circle = 10, n.prev = 3)
	}

	sparsity <- sum(X == 0) / prod(dim(X))
	cat(sprintf('[%s] sparsity: %.3f%%\n', Sys.time(), sparsity * 100))

	if (missing(init)){
		if (sparsity < 0.6 && ncol(CT) <= 5)
			init <- list(method = 'all', update.beta = TRUE)
		else if (sparsity < 0.6 && ncol(CT) > 5)
			init <- list(method = 'forward', update.beta = TRUE)
		else if (sparsity >= 0.6 && ncol(CT) <= 5)
			init <- list(method = 'backward', update.beta = FALSE)
		else if (sparsity >= 0.6 && ncol(CT) > 5)
			init <- list(method = 'forward', update.beta = FALSE)
	}
	cat(sprintf('[%s] gtm initialization method: %s\n', Sys.time(), init$method))
	
	K <- ls$K
	cat(sprintf('[%s] number of input rows(N): %d\n', Sys.time(), N))
	cat(sprintf('[%s] number of input colnums(M): %d\n', Sys.time(), M))
	cat(sprintf('[%s] number of time points(ncol(time.table)): %d\n', Sys.time(), ncol(CT)))
	cat(sprintf('[%s] number of metagenes(K): %d\n', Sys.time(), K))
	cat(sprintf('[%s] optimization method: %s\n', Sys.time(), control[['optimization.method']]))
	if (control[['optimization.method']] == 'stochastic'){
		cat(sprintf('[%s] batch.size: %d\n', Sys.time(), control[['batch.size']]))
	}
	cat(sprintf('[%s] number of cores: %d\n', Sys.time(), control[['mc.cores']]))

	cat(sprintf('[%s] initializing metagene coefficients (max.size.dist.matrix=%d):\n', Sys.time(), control$max.size.dist.matrix))
	V <- fastmds(scale(log(X + 1)), K = K, scale = TRUE, control = control)

	cat(sprintf('[%s] initializing metagene basis:\n', Sys.time()))
	a <- rep(0, M)	# initialize the cell effect
	CC <- sparseMatrix(i = 1:M, j = pam(t(V), K, cluster.only = TRUE), dims = c(M, K))	# # cells ~ clusters
	U <- as.matrix(X %*% CC %*% diag(1 / Matrix::colSums(CC)))	# initialize E[U]
	mf <- ccm(X, U = U, V = V, a = a, control = control)
	U <- mf$U	# the metagene basis
	V <- mf$V	
	a <- mf$a # the cell effect

	cat(sprintf('[%s] initializing prototype landscape (method=%s):\n', Sys.time(), init$method))
	mf <- gtm(V = V, landscape = ls, CT = CT, method = init$method, update.beta = init$update.beta, control = control)

	mem <- max.col(mf$Z)
	mem[apply(mf$Z, 1, max) < mf$landscape$assignment.threshold] <- NA
	mf <- tcm.core(X = X, U = U, V = V, a = a, landscape = mf$landscape, CT = CT, mem = mem, control = control)

	mf

} # end of tcm


#' softmax
#' 
#' Computing y[i, j] = exp(x[i, j]) / sum(exp(x[, j]))
#'
softmax <- function(x){
	x.max <- apply(x, 2, max)
	y <- log(Matrix::rowSums(exp(t(x) - x.max)) + .Machine$double.eps) + x.max	# log(sum(exp(x[, j])))
	y <- t(exp(t(x) - y))
	y
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
tcm.core <- function(X, U = NULL, V = NULL, a = NULL, landscape = NULL, CT = NULL, mem = NULL, c0 = 0.001, d0 = 0.001, control = NULL){

	eps <- .Machine$double.eps
	control <- set.control(control)

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- ncol(U)	# number of metagenes
	H <- landscape[['H.prototype']]	# number of prototypes

	if (!is.null(landscape)){
		Theta <- as.matrix(landscape[['Theta.free']] %*% landscape[['S']])	# all prototypes
		if (landscape$type == 'temporal.convolving')
			W <- CT %*% landscape[['TMC']]	# a binary cell ~ prototype assignment matrix
		else
			W <- matrix(TRUE, M, landscape[['H.prototype']])
	}else
		stop('landscape must be specified')

	# pre-assigned cells can be only mapped to dedicated prototypes
	if (!is.null(mem) && any(!is.na(mem))){
		i <- !is.na(mem)
		W[i, ] <- FALSE
		W[cbind(which(i), mem[i])] <- TRUE
	}

	if (is.null(a))
		a <- rep(0, M)	# the cell effect

	# splitting the samples into batches
	s <- split.data(M, control = control)

	# initialize the local variables
	localvar <- mclapply(1:s[['n.batch']], function(b){
		m <- s[['groups']] == b
		V.var <- matrix(1, nrow = K, ncol = s[['size']][b])	# Var[V]
		V.exp <- exp(V[, m] + 1/2 * V.var)	# E[exp(V)]
		list(a = a[m], V = V[, m], V.exp = V.exp, Z = matrix(1 / H, s[['size']][b], H))
	}, mc.cores = control[['mc.cores']])

	iter <- 1
	optval <- NULL
	while (iter <= control[['max.iter']]){

		Theta.free.p <- landscape[['Theta.free']]

		# updating all local variables: E[V], E[exp(V)], a and E[Z]
		localvar <- mclapply(1:s[['n.batch']], function(b){
			m <- s[['groups']] == b
			update.localvar(
				X = X[, m, drop = FALSE],
				U = U,
				V = localvar[[b]]$V, 
				V.exp = localvar[[b]]$V.exp,
				W = W[m, ],
				landscape = landscape,
				model = 'tcm', control = control  
			)
		}, mc.cores = control[['mc.cores']])

		rho <- learning.rate(iter, control)
		m <- sample.data(s, control)

		# updating free prototypes Theta
		SZS <- landscape[['S']] %*% Diagonal(x = landscape[['beta']] * Reduce('+', mclapply(1:s[['n.batch']], function(b) colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]), mc.cores = control[['mc.cores']]))) %*% t(landscape[['S']])
		SZSL.inv <- chol2inv(chol(as.matrix(SZS + landscape[['L']])))
		landscape[['Theta.free']] <- (1 - rho) * landscape[['Theta.free']] + rho * landscape[['beta']] * Reduce('+', mclapply(1:s[['n.batch']], function(b){
			localvar[[b]]$V[, m[[b]], drop = FALSE] %*% localvar[[b]]$Z[m[[b]], , drop = FALSE] %*% t(landscape[['S']]) %*% SZSL.inv
		}, mc.cores = control[['mc.cores']]))

		U <- (1 - rho) * U + rho * Reduce('+', mclapply(1:s[['n.batch']], function(b){
			c <- U * (X[, which(s[['groups']] == b)[m[[b]]]] / (U %*% localvar[[b]]$V.exp[, m[[b]], drop = FALSE] + eps)) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE]) + c0
			d <- matrix(1, N, length(m[[b]])) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE] %*% diag(exp(localvar[[b]]$a[m[[b]]]))) + d0
			U.log <- digamma(c) - log(d + eps)	# E[log(U)]
			exp(U.log)
		}, mc.cores = control[['mc.cores']])) / s[['n.batch']]

		landscape[['beta']] <- (1 - rho) * landscape[['beta']] + rho * Reduce('+', mclapply(1:s[['n.batch']], function(b){
			(K * sum(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + 0.001) / (sum(localvar[[b]]$Z[m[[b]], , drop = FALSE] * rdist(t(localvar[[b]]$V[, m[[b]], drop = FALSE]), t(Theta))^2) + 0.001)
		}, mc.cores = control[['mc.cores']])) / s[['n.batch']]

		# updating kappa/pi.log
		kappa <- Reduce('+', mclapply(1:s[['n.batch']], function(b){
			colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + landscape[['alpha0']]
		}, mc.cores = control[['mc.cores']])) 
		landscape[['pi.log']] <- (1 - rho) * landscape[['pi.log']] + rho * (digamma(kappa) - digamma(sum(kappa))) # E[log(pi)]

		if (iter == 1 || iter %% 10 == 0){
			J <- norm(landscape[['Theta.free']] - Theta.free.p, '2')
			n.cluster <- length(unique(unlist(mclapply(1:s[['n.batch']], function(b) max.col(localvar[[b]]$Z), mc.cores = control[['mc.cores']]))))
			cat(sprintf('[%s] tcm | iter=%5.d | M=%5.d | n.batch=%3.d | dJ=%7.3e | # non-empty prototypes=%4.d | beta=%7.3f\n', Sys.time(), iter, M, s[['n.batch']], J, n.cluster, landscape[['beta']]))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), n.cluster = n.cluster, J = J))
			if (J < 0.01)
				break
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	a <- rep(0, M)
	V <- matrix(0, K, M)
	Z <- matrix(0, M, H)
	for (b in 1:s[['n.batch']]){
		V[, s[['groups']] == b] <- localvar[[b]]$V
		Z[s[['groups']] == b, ] <- localvar[[b]]$Z
		a[s[['groups']] == b] <- localvar[[b]]$a
	}

	cat(sprintf('[%s] tcm | updating active prototypes\n', Sys.time()))
	landscape$is.active <- 1:landscape$H.prototype %in% max.col(Z)

	cat(sprintf('[%s] tcm | updating landscape graph\n', Sys.time()))
	landscape$graph <- as.igraph(landscape)

	structure(list(H = H, N = N, M = M, K = K, U = U, V = V, a = a, Z = Z, CT = CT, landscape = landscape, optval = optval), class = 'tcm')

} # end of tcm.core



#' plot.tcm
#' 
#' Plot the TCM results
#' 
#' @param x a tcm object
#' @param ... Further graphical parameters may also be supplied as arguments
#'
#' @export
#'
#' @author Wuming Gong, \email{gongx030@umn.edu}
#'
#' @examples
#' library(tcm)
#' set.seed(122)
#' sim <- sim.rnaseq.ts(N = 2000, M = 500, n.lineage = 5, n.time.points = 5)
#' X <- assays(sim)$count
#' time.table <- colData(sim)$time.table
#' mf <- tcm(X, time.table = time.table)
#' bg.cell <- rainbow(5)[colData(sim)$lineage]
#' dev.new(height = 10, width = 12)
#' par(mar = c(5, 5, 5, 15))
#' plot(mf, pch = 21, bg = bg.cell, cex = 1.5)
#' 
plot.tcm <- function(x, ...){

	param <- list(...)

	if (x$landscape$type %in% c('temporal.convolving', 'plate')){

		H.prototype <- x$landscape$H.prototype	# number of prototypes
		M <- ncol(x$V)	# number of cells
		mem <- max.col(x$Z)	# cluster membership
		coord <- x$landscape$Y.prototype
		csize <- table(factor(mem, 1:H.prototype))	# size of cells that were assigned to each prototype
		ncls <- sum(csize > 0)	# number of prototypes that have at least one assigned cells
		centers <- which(csize > 0)	# index of non-empty prototype
		CT <- x$CT
		Yp <- x$landscape$Yp.prototype
		Y <- x$landscape$Y.prototype
		n.prototype <- x$landscape$n.prototype

		lim <- c(-max(Yp[, 'r']), max(Yp[, 'r']))	# limit of the 2D space
		plot(NA, xlim = lim * 1.1, ylim = lim * 1.1, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', asp = 1, main = param[['main']])

		rs <- unique(Yp[, 'r'])
		col.time <- colorpanel(length(rs), low = 'gray100', mid = 'gray70', high = 'gray40')
		for (i in length(rs):1){
			draw.circle(0, 0, rs[i], nv = 100, border = NA, col = col.time[i])
		}
		
		if (x[['landscape']][['type']] == 'temporal.convolving'){
			for (t in 1:ncol(x[['CT']]))
				draw.circle(0, 0, t, nv = 100, border = 'black', lty = 2, lwd = 2)
		}

		Y.cell2 <- t(jitter2(t(coord[, max.col(x$Z)]), amount = 0.125 * max(lim)))	# jittered coordinate of cells
		points(Y.cell2[1, ], Y.cell2[2, ], ...)
	}

} # end of plot.tcm


#' Add the labels of active prototypes
#' 
centers <- function(mf, ...){
	coord <- mf$landscape$Y.prototype
	c <- which(mf$landscape$is.active)
	text(coord[1, c], coord[2, c], c, ...)
} # end of centers


#' Adding an Enclidean noise, x must be in a two dimensional matrix
#'
jitter2 <- function(x, amount = 0.1){
	phi <- runif(nrow(x), min = 0, max = 2 * pi)
	r <- runif(nrow(x), min = 0, max = amount)
	x[, 1] <- x[, 1] + r * cos(phi)
	x[, 2] <- x[, 2] + r * sin(phi)
	x
} # end of jitter2


#' num2color
#'
#' Scale a continuous vector into color representation
#' @param x a numeric vector
#' @param cp the color spectrum (default: colorpanel(100, low = 'blue', mid = 'black', high = 'yellow'))
#'
#' @return a character vector of HTML-style RGB colors 
#'
#' @export
#'
#' @author Wuming Gong, \email{gongx030@umn.edu}
#'
num2color <- function(x, cp = colorpanel(100, low = 'blue', mid = 'black', high = 'yellow')) cp[round((x - min(x)) / (max(x) - min(x)) * (length(cp) - 1)) + 1]


#' learning rate for stochastic variational inference
#'
#' @param t current iteration step
#' @param decay.rate down-weights for early iterations (default: 5)
#' @param forgetting.rate controls how quickly old information is forgotten (default: 0.75)
#'
learning.rate <- function(t, control){
	rho <- switch(control[['optimization.method']],
		'stochastic' = learning.rate.stochastic(t, control[['decay.rate']], control[['forgetting.rate']]),
		'batch' = 1
	)
	rho
}

learning.rate.stochastic <- function(t, decay.rate, forgetting.rate) ifelse(t == 0, 1, (t + decay.rate)^(-forgetting.rate))


#' updating local variable E[V], E[exp(V)], a and E[Z]
#' @param X a N by M expression matrix
#' @param U a N by K metagene basis matrix
#' @param V the initial values of a K by M metagene coefficient matrix
update.localvar <- function(X, U, V = NULL, V.exp = NULL, W = NULL, landscape = NULL, model = NA, control = NULL, ...){

	eps <- .Machine$double.eps
	param <- list(...)

	max.iter <- switch(control[['optimization.method']], 'batch' = 50, 'stochastic' = 1)

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
			Theta <- as.matrix(landscape[['Theta.free']] %*% landscape[['S']])	# all prototypes
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
		Theta <- landscape$Theta.free %*% landscape$S
		Z <- as.matrix(-1/2  * landscape$beta * rdist(t(V), t(Theta))^2 - 1/2 * K * log(2 * pi) + 1/2 * K * log(landscape$beta) + matrix(landscape$pi.log, M, H, byrow = TRUE))
		if (!is.null(W))
			Z[!W] <- -Inf
		Z <- t(softmax(t(Z)))# E[Z]
		list(V = V, Z = Z)
	}
} # end of update.localvar


#' landscape
#' 
#' Define the prototype landscape for simulating and visualizing scRNA-seq data 
#' 
#' @param type type of the prototype landscape (default: temporal.convolving)
#' 
#' @param n.prototype prototypes per layer
#' @param n.circle the number of layer 
#' @param n.prev number of convolving layers per time point
#' 
#' @return a landscape object
#'
#' @export
#'
#' @author Wuming Gong, \email{gongx030@umn.edu}
#'
#' @examples
#' # a temporal convolving landscape with the following parameters:
#' #  number of metagenes (K): 15
#' #  number of prototypes per layer (n.prototype): 15
#' #  number of layers per time point (n.circle): 10
#' #  number of convolving layers per time point (n.prev): 3
#' ls <- landscape(type = 'temporal.convolving', K = 15, n.prototype = 15, n.circle = 10, n.prev = 3)
#'
#' # a plate landscape with following parameters:
#' #  number of metagenes (K): 15
#' #  number of prototypes per layer (n.prototype): 10
#' #  number of layers (n.circle): 100
#' ls <- landscape(type = 'plate', K = 15, n.prototype = 10, n.circle = 100)
#'
landscape <- function(type = 'temporal.convolving', K = 15, ...){

	param <- list(...)

	if (is.null(param[['alpha00']]))
		alpha00 <- 0.001
	else
		alpha00 <- param[['alpha00']]

	if (type == 'temporal.convolving'){

		time.points <- param[['time.points']]

		if (is.null(time.points) || time.points < 2)
			stop('at least two time points must be specified for landscape tcm')

		if (is.null(param$n.prototype)){
			param$n.prototype <- 15
			cat(sprintf('[%s] number of prototypes per circle: %d\n', Sys.time(), param$n.prototype))
		}

		if (is.null(param$n.circle)){
			param$n.circle <- 10
			cat(sprintf('[%s] number of circles per time point: %d\n', Sys.time(), param$n.circle))
		}

		if (is.null(param$n.prev)){
			param$n.prev <- min(3, param$n.circle - 1)
			cat(sprintf('[%s] number of circle(s) mapped from the previous time point: %d\n', Sys.time(), param$n.prev))
		}

		n.prototype <- param$n.prototype
		n.circle <- param$n.circle
		n.prev <- param$n.prev

		if (is.null(param$lambda))
			lambda <- 1
		else
			lambda <- param$lambda

		H <- n.circle * n.prototype # number of prototypes per time point
		H.prototype <- H * time.points	# total number of prototypes
		alpha0 <- rep(alpha00, H.prototype)

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

		for (tp in 2:time.points){
			hc <- TMC[tp, ] 	# prototypes from current time point
			hp <- TMC[tp - 1, ] 	# prototypes from previous time point
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
		alpha0 <- rep(alpha00, H.prototype)
		pi.log <- digamma(alpha0) - digamma(sum(alpha0))	# E[log(pi)]

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
		Theta.free <- mvrnorm(K, rep(0, H.prototype), Sigma = L.inv)
		alpha0 <- rep(alpha00, H.prototype)
		pi.log <- digamma(alpha0) - digamma(sum(alpha0))	# E[log(pi)]

	}else if  (type == 'ladder'){

		if (is.null(param[['lambda']]))
			lambda <- 1
		else
			lambda <- param[['lambda']]

		time.points <- param$time.points
		if (is.null(time.points) || time.points < 2)
			stop('at least two time points must be specified for ladder landscape')

		n.prototype <- param$n.prototype
		n.circle <- param$n.circle
		H.free <- H.prototype <- H <- time.points * n.circle * n.prototype # number of prototypes 
		Y.prototype <- t(as.matrix(expand.grid(x = seq(-1, 1, length.out = n.prototype), y = seq(-1, 1, length.out = n.circle * time.points))))
		Yp.prototype <- NULL
		L.inv <- rdist(t(Y.prototype))	# pairwise euclidean distance between prototypes
		L.inv <- exp(-L.inv^2 / (2 * lambda^2))	# squared exponential covariance function
		L.inv <- as.matrix(nearPD(L.inv)$mat)	# convert to a nearest PSD matrix
		L <- chol2inv(chol(L.inv))	# inverse the corrected covariance matrix by Cholesky decomposition
		S <- Diagonal(n = H.prototype)
		is.free <- rep(TRUE, H.prototype)
		n.prev <- NULL
		TMC <- sparseMatrix(i = rep(1:time.points, each = n.circle * n.prototype), j = 1:H.prototype, dims = c(time.points, H.prototype))	# time ~ prototypes
		MC <- NULL
		Theta.free <- mvrnorm(K, rep(0, H.prototype), Sigma = L.inv)
		alpha0 <- rep(alpha00, H.prototype)
		pi.log <- digamma(alpha0) - digamma(sum(alpha0))	# E[log(pi)]

	}else if (type == 'temporal.plate'){

		if (is.null(param$lambda))
			lambda <- 1
		else
			lambda <- param$lambda

		time.points <- param$time.points
		if (is.null(time.points) || time.points < 2)
			stop('at least two time points must be specified for temporal plate landscape')

		n.prototype <- param$n.prototype
		n.circle <- param$n.circle
		H.free <- H.prototype <- H <- time.points * n.circle * n.prototype # number of prototypes 
		Yp.prototype <- as.matrix(expand.grid(angle = seq(0, 2 * pi, length.out = n.prototype + 1)[-1], r = seq(0, 1, length.out = n.circle * time.points + 1)[-1]))	# polar coord
		Y.prototype <- rbind(x = Yp.prototype[, 'r'] * cos(Yp.prototype[, 'angle']), y = Yp.prototype[, 'r'] * sin(Yp.prototype[, 'angle']))	# x-y coord
		L.inv <- rdist(t(Y.prototype))	# pairwise euclidean distance between prototypes
		L.inv <- exp(-L.inv^2 / (2 * lambda^2))	# squared exponential covariance function
		L.inv <- as.matrix(nearPD(L.inv)$mat)	# convert to a nearest PSD matrix
		L <- chol2inv(chol(L.inv))	# inverse the corrected covariance matrix by Cholesky decomposition
		S <- Diagonal(n = H.prototype)
		is.free <- rep(TRUE, H.prototype)
		n.prev <- NULL
		TMC <- sparseMatrix(i = rep(1:time.points, each = n.circle * n.prototype), j = 1:H.prototype, dims = c(time.points, H.prototype))	# time ~ prototypes
		MC <- NULL
		Theta.free <- mvrnorm(K, rep(0, H.prototype), Sigma = L.inv)
		alpha0 <- rep(alpha00, H.prototype)
		pi.log <- digamma(alpha0) - digamma(sum(alpha0))	# E[log(pi)]

	}else
		stop(sprintf('unknown landscape type: %s', type))

	landscape <- structure(list(
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
		alpha00 = alpha00,
		alpha0 = alpha0,
		is.free = is.free,
		pi.log = pi.log,
		beta = 1,
		is.active = rep(TRUE, H.prototype),
		is.mappable = rep(TRUE, H.prototype),
		assignment.threshold = 0.95
	), class = 'landscape')

	landscape$graph <- as.igraph(landscape)
	landscape

} # end of landscape 


#' Get an igraph object of neighboring graph of the prototype landscape where each vertex represent a prototype
#' and the distance between neighboring vertex is the Euclidean distance of their metagene coefficients
#'
#' @param x a tcm object
#'
as.igraph.landscape <- function(x, ...){

	if (x$type %in% c('temporal.convolving', 'plate', 'temporal.plate'))
		Yp <- transform(x$Yp.prototype, angle = factor(angle / (2 * pi)), r = factor(r))
	else if (x$type %in% c('ladder'))
		Yp <- data.frame(angle = factor(x$Y.prototype[1, ]), r = factor(x$Y.prototype[2, ]))

	rs <- levels(Yp[, 2])
	angles <- levels(Yp[, 1])
	a <- rbind( 
		do.call('rbind', lapply(rs, function(r){
			v <- which(Yp[, 2] == r)
			data.frame(from = v[1:length(v)], to = v[c(2:length(v), 1)])
		})),
		do.call('rbind', lapply(angles, function(angle){
			v <- which(Yp[, 1] == angle)
			data.frame(from = v[1:(length(v) - 1)], to = v[2:length(v)])
		}))
	)
		
	Theta <- as.matrix(x[['Theta.free']] %*% x[['S']])	# all prototypes
	w <- sqrt(Matrix::colSums((Theta[, a[, 1]] - Theta[, a[, 2]])^2))	# distance between neighboring prototypes
	a <- cbind(a, weight = w)
	A <- sparseMatrix(i = a[, 'from'], j = a[, 'to'], x = a[, 'weight'], dims = c(x$H.prototype, x$H.prototype))
	A <- A + t(A)
	A <- graph.adjacency(A, weighted = TRUE, mode = 'undirected')
	A
} # end of as.igraph.landscape


#' Correlated Cell Model
#'
ccm <- function(X, U = NULL, V = NULL, a = NULL, c0 = 0.001, d0 = 0.001, control = NULL){

	eps <- .Machine$double.eps
	control <- set.control(control)

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- nrow(V)	# number of metagenes

	# splitting the samples into batches
	s <- split.data(M, control = control)

	# initialize the local variables
	localvar <- mclapply(1:s[['n.batch']], function(b){
		m <- s[['groups']] == b
		V.var <- matrix(1, nrow = K, ncol = s[['size']][b])	# Var[V]
		V.exp <- exp(V[, m] + 1/2 * V.var)	# E[exp(V)]
		list(a = a[m], V = V[, m], V.exp = V.exp)
	}, mc.cores = control[['mc.cores']])

	iter <- 1
	optval <- NULL
	while (iter <= control[['ccm.max.iter']]){

		Up <- U

		# updating all local variables: E[V], E[exp(V)] and a
		localvar <- mclapply(1:s[['n.batch']], function(b){
			m <- s[['groups']] == b
			update.localvar(
				X = X[, m, drop = FALSE],
				U = U,
				V = localvar[[b]]$V, 
				V.exp = localvar[[b]]$V.exp,
				model = 'ccm', control = control
			)
		}, mc.cores = control[['mc.cores']])

		rho <- learning.rate(iter, control)
		m <- sample.data(s, control)

		U <- (1 - rho) * U + rho * Reduce('+', mclapply(1:s[['n.batch']], function(b){
			c <- U * (X[, which(s[['groups']] == b)[m[[b]]]] / (U %*% localvar[[b]]$V.exp[, m[[b]], drop = FALSE] + eps)) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE]) + c0
			d <- matrix(1, N, length(m[[b]])) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE] %*% diag(exp(localvar[[b]]$a[m[[b]]]))) + d0
			U.log <- digamma(c) - log(d + eps)	# E[log(U)]
			exp(U.log)
		}, mc.cores = control[['mc.cores']])) / s[['n.batch']]

		if (iter == 1 || iter %% 10 == 0){

			J.X <- sum(unlist(mclapply(1:s[['n.batch']], function(b) sum(ldpois(X[, which(s[['groups']] == b)], U %*% localvar[[b]]$V.exp %*% diag(exp(localvar[[b]]$a)))), mc.cores = control[['mc.cores']])))
			J.V <- -sum(unlist(mclapply(1:s[['n.batch']], function(b) norm(localvar[[b]]$V, '2')^2, mc.cores = control$mc.cores))) / (K * M)
			cat(sprintf('[%s] ccm | iter=%5.d | M=%5.d | n.batch=%3.d | J(X)=%7.3e | J(V)=%7.3e | J=%7.3e\n', Sys.time(), iter, M, s[['n.batch']], J.X, J.V, J.X + J.V))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), J.X = J.X, J.V = J.V, J = J.X + J.V))
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	a <- rep(0, M)
	V <- matrix(0, K, M)
	for (b in 1:s[['n.batch']]){
		V[, s[['groups']] == b] <- localvar[[b]]$V
		a[s[['groups']] == b] <- localvar[[b]]$a
	}
	list(N = N, M = M, K = K, U = U, V = V, a = a, optval = optval)

} # end of ccm


#' The generative topographic model
#'
#' @param V a K by M metagene coefficient of cells, where K is the number of metagenes and M is the number of cells
gtm.core <- function(V, landscape = NULL, CT = NULL, mem = NULL, update.beta = FALSE, control = NULL){

	eps <- .Machine$double.eps
	M <- ncol(V)	# number of cells
	K <- nrow(V)	# number of metagenes
	control <- set.control(control)

	if (!is.null(landscape)){
		Theta <- as.matrix(landscape$Theta.free %*% landscape$S)	# all prototypes
		if (landscape[['type']] %in% c('temporal.convolving', 'temporal.plate')){
			if (is.null(CT))
				stop('CT cannot be NULL') 
			W <- CT %*% landscape$TMC	# a binary cell ~ prototype assignment matrix
		}else if (landscape[['type']] %in% c('som', 'plate')){
			W <- as(matrix(TRUE, M, landscape$H.prototype), 'ngCMatrix')
		}else if (landscape$type %in% c('ladder')){
			W <- CT %*% landscape$TMC	# a binary cell ~ prototype assignment matrix
		}else
			stop(sprintf('unknown landscape: %s', landscape$type))
	}else
		stop('landscape must be specified')

	# turn off the non-mappable prototypes
	if (any(!landscape$is.mappable))
		W[, !landscape$is.mappable] <- FALSE
	
	# pre-assigned cells can be only mapped to dedicated prototypes
	if (!is.null(mem) && any(!is.na(mem))){
		i <- !is.na(mem)
		W[i, ] <- FALSE
		W[cbind(which(i), mem[i])] <- TRUE
	}
	
	# splitting the samples into batches
	s <- split.data(M, control = control)


	iter <- 1
	optval <- NULL
	while (iter <= control[['max.iter']]){

		Theta.free.p <- landscape[['Theta.free']]

		# updating all local variables: E[Z]
		localvar <- mclapply(1:s[['n.batch']], function(b){
			m <- s[['groups']] == b
			update.localvar(
				V = V[, m, drop = FALSE],
				W = W[m, , drop = FALSE],
				landscape = landscape,
				model = 'gtm', control = control
			)
		}, mc.cores = control[['mc.cores']])

		rho <- learning.rate(iter, control)
		m <- sample.data(s, control)

		# updating free prototypes Theta
		SZS <- landscape$S %*% Diagonal(x = landscape$beta * Reduce('+', mclapply(1:s$n.batch, function(b) colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]), mc.cores = control$mc.cores))) %*% t(landscape$S)
		SZSL.inv <- chol2inv(chol(as.matrix(SZS + landscape$L)))
		landscape$Theta.free <- (1 - rho) * landscape$Theta.free + rho * landscape$beta * Reduce('+', mclapply(1:s$n.batch, function(b){
			Vb <- localvar[[b]]$V[, m[[b]], drop = FALSE]																																																					 
			Zb <- localvar[[b]]$Z[m[[b]], , drop = FALSE]
			Vb %*% Zb %*% t(landscape$S) %*% SZSL.inv
		}, mc.cores = control$mc.cores))
		Theta <- as.matrix(landscape$Theta.free %*% landscape$S)	# all prototypes

		if (update.beta){
			landscape$beta <- (1 - rho) * landscape[['beta']] + rho * Reduce('+', mclapply(1:s[['n.batch']], function(b){
				(K * sum(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + 0.001) / (sum(localvar[[b]]$Z[m[[b]], , drop = FALSE] * rdist(t(localvar[[b]]$V[, m[[b]], drop = FALSE]), t(Theta))^2) + 0.001)
			}, mc.cores = control[['mc.cores']])) / s[['n.batch']]
		}

		# updating pi.log
		kappa <- Reduce('+', mclapply(1:s[['n.batch']], function(b){
			colSums(localvar[[b]]$Z[m[[b]], , drop = FALSE]) + landscape$alpha0
		}, mc.cores = control[['mc.cores']])) 
		landscape[['pi.log']] <- (1 - rho) * landscape[['pi.log']] + rho * (digamma(kappa) - digamma(sum(kappa))) # E[log(pi)]
		
		if (iter == 1 || iter %% 10 == 0){
			J <- norm(landscape[['Theta.free']] - Theta.free.p, '2')
			n.cluster <- length(unique(unlist(mclapply(1:s[['n.batch']], function(b) max.col(localvar[[b]]$Z), mc.cores = control$mc.cores))))
			Mc <- sum(unlist(mclapply(1:s$n.batch, function(b) apply(localvar[[b]]$Z, 1, max) > landscape$assignment.threshold, mc.cores = control$mc.cores)))
			cat(sprintf('[%s] gtm | iter=%5.d | M=%d | n.batch=%d | Mc=%d | dJ=%7.3e | # non-empty prototypes=%4.d | beta=%7.3f\n', Sys.time(), iter, M, s$n.batch, Mc, J, n.cluster, landscape[['beta']]))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), n.cluster = n.cluster, beta = landscape[['beta']], J = J))
			if (J < 0.01)
				break
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	Z <- matrix(0, M, landscape[['H.prototype']], dimnames = list(colnames(V), NULL))
	for (b in 1:s[['n.batch']]){
		Z[s[['groups']] == b, ] <- localvar[[b]]$Z
	}

	cat(sprintf('[%s] gtm | updating active prototypes\n', Sys.time()))
	landscape$is.active <- 1:landscape$H.prototype %in% max.col(Z)

	cat(sprintf('[%s] gtm | updating landscape graph\n', Sys.time()))
	landscape$graph <- as.igraph(landscape)

	structure(list(M = M, K = K, V = V, Z = Z, CT = CT, optval = optval, landscape = landscape), class = 'gtm')

} # end of gtm


#' log density of Poission distribution
#' (faster than dpois(..., log = TRUE) due to the matrix operations)
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
		cat(sprintf('number of time points: %d\n', nrow(x[['TMC']])))
		cat(sprintf('number of circle per time point: %d\n', x[['n.circle']]))
		cat(sprintf('number of prototypes per circle: %d\n', x[['n.prototype']]))
		cat(sprintf('number of prototypes for each time point: %d\n', x[['H']]))
		cat(sprintf('number of circle(s) mapped from the previous time point: %d\n', x[['n.prev']]))
	}else if (x[['type']] == 'temporal.plate'){
		cat(sprintf('number of time points: %d\n', nrow(x[['TMC']])))
		cat(sprintf('number of circle per time point: %d\n', x[['n.circle']]))
		cat(sprintf('number of prototypes per circle: %d\n', x[['n.prototype']]))
		cat(sprintf('number of prototypes for each time point: %d\n', x[['H']]))
	}else if (x[['type']] == 'plate'){
		cat(sprintf('number of circles: %d\n', x[['n.circle']]))
		cat(sprintf('number of prototypes per circle: %d\n', x[['n.prototype']]))
	}
} # end of print.landscape


#' Correlated Cell Regression
#'
ccr <- function(X, U = NULL, V = NULL, a = NULL, c0 = 0.001, d0 = 0.001, control = NULL){

	eps <- .Machine$double.eps
	control <- set.control(control)

	N <- nrow(X)	# number of genes
	M <- ncol(X)	# number of cells
	K <- nrow(V)	# number of metagenes

	# splitting the samples into batches
	s <- split.data(M, control = control)

	# initialize the local variables
	localvar <- mclapply(1:s[['n.batch']], function(b) list(a = a[s[['groups']] == b]), mc.cores = control[['mc.cores']])

	if (is.null(U))
		stop('U cannot be NULL')

	iter <- 1
	optval <- NULL
	while (iter <= control[['max.iter']]){

		# updating all local variables: E[V], E[exp(V)] and a
		localvar <- mclapply(1:s[['n.batch']], function(b){
			m <- s[['groups']] == b
			update.localvar(
				X = X[, m, drop = FALSE],
				U = U,
				V = V[, m, drop = FALSE],
				model = 'ccr', control = control
			)
		}, mc.cores = control[['mc.cores']])

		rho <- learning.rate(iter, control)
		m <- sample.data(s, control)

		U <- (1 - rho) * U + rho * Reduce('+', mclapply(1:s[['n.batch']], function(b){
			c <- U * (X[, which(s[['groups']] == b)[m[[b]]]] / (U %*% localvar[[b]]$V.exp[, m[[b]], drop = FALSE] + eps)) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE]) + c0
			d <- matrix(1, N, length(m[[b]])) %*% t(localvar[[b]]$V.exp[, m[[b]], drop = FALSE] %*% diag(exp(localvar[[b]]$a[m[[b]]]))) + d0
			U.log <- digamma(c) - log(d + eps)	# E[log(U)]
			exp(U.log)
		}, mc.cores = control[['mc.cores']])) / s[['n.batch']]

		if (iter == 1 || iter %% 10 == 0){
			J <- sum(unlist(mclapply(1:s[['n.batch']], function(b) sum(ldpois(X[, which(s[['groups']] == b)], U %*% localvar[[b]]$V.exp %*% diag(exp(localvar[[b]]$a)))), mc.cores = control[['mc.cores']])))
			cat(sprintf('[%s] ccr | iter=%5.d | M=%5.d | n.batch=%3.d | J=%7.3e\n', Sys.time(), iter, M, s[['n.batch']], J))
			optval <- rbind(optval, data.frame(iter = iter, time = Sys.time(), J = J))
		}
		iter <- iter + 1
	}

	# aggregating the results from different batches
	a <- rep(0, M)
	for (b in 1:s[['n.batch']]){
		a[s[['groups']] == b] <- localvar[[b]]$a
	}
	list(N = N, M = M, K = K, U = U, a = a, optval = optval, control = control)

} # end of ccr


sample.data <- function(s, control){
	if (control[['optimization.method']] == 'stochastic'){
		# SVI: sampling a subset of samples(cells) for updating global variables
		m <- lapply(1:s[['n.batch']], function(b) sample(1:s[['size']][b], min(max(1, round(control[['batch.size']] / s[['n.batch']])), s[['size']][b])))
	}else{
		# batch optimization: use all the samples
		m <- lapply(1:s[['n.batch']], function(b) 1:s[['size']][b])
	}
	m
} # end of smaple.data


#' Splitting data into different batches
#' 
split.data <- function(M, f = NULL, control){

	if (is.null(f)){	# if there no pre-defined groups
		if (M <= control$mc.cores * control$max.size.per.batch){	
			n.batch <- control$mc.cores
		}else{
			n.batch <- ceiling(M / control$max.size.per.batch)
		}
		groups <- sample(1:n.batch, M, replace = TRUE)
	}else{
		sp <- split(1:M, list(factor(f)))	
		groups <- lapply(1:length(sp), function(i){
			Mi <- length(sp[[i]])
			if (Mi > control$max.size.per.batch){
				n.batch.i <- ceiling(Mi / control$max.size.per.batch)
				sample(1:n.batch.i, Mi, replace = TRUE)
			}else
				rep(1, Mi)
		})
		groups <- as.numeric(factor(sprintf('%d_%d', rep(1:length(sp), sapply(sp, length)), unlist(groups))))
		n.batch <- max(groups)
	}
	size <- table(factor(groups, 1:n.batch))	# number of samples per core

	# splitting the samples into batches
	list(size = size, n.batch = n.batch, groups = groups)

} # end of split.data

#' Setting the default values for the arguments for optimizing TCM
#'
#' @param max.size.dist.matrix The maximum size of the distance matrix (for fastmds)

set.control <- function(control){
	control.default <- list(
		optimization.method = 'batch', 
		batch.size = 2000, 
		ccm.max.iter = 10,
		max.iter = 50, 
		mc.cores = 1, 
		decay.rate = 1, 
		forgetting.rate = 0.75, 
		max.size.per.batch = 2000,
		max.size.dist.matrix = 2000,
		fine.tune = FALSE
	)
	for (i in names(control.default))
		if (!is.null(control[[i]]))
			control.default[[i]] <- control[[i]]
	control.default	
} # end of set.control


#' Plot the paths on the prototype landscape
#' @param x a tcm object
#' @param ...  Arguments to be passed to arrow(s)
#'
add.paths.landscape <- function(x, ...){

	param <- list(...)

	if (is.null(x$paths)){
		stop('landscape$path is NULL')
	}else{
		Y <- x$Y.prototype
		g <- do.call('rbind', lapply(1:length(x$paths), function(i){
			vs <- x$paths[[i]]
			cbind(from = vs[1:(length(vs) - 1)], to = vs[2:length(vs)])
		}))
		g <- unique(g)
		segments(Y[1, g[, 'from']], Y[2, g[, 'from']], Y[1, g[, 'to']], Y[2, g[, 'to']], ...)
	}

} # end of add.paths.landscape

add.paths.tcm <- function(x, ...) add.paths(x$landscape, ...)


#' Find the trajectory
trajectory <- function(x, ...){

	Y <- x$landscape$Y.prototype	# 2D coordinates of prototypes on the landscape

	if (x$landscape$type %in% c('temporal.convolving', 'plate')){

		mem <- max.col(x$Z)
		centers <- which(x$landscape$is.active)
		is.center <- 1:x$landscape$H.prototype %in% centers
		Theta <- as.matrix(x$landscape$Theta.free %*% x$landscape$S)	# metagene coefficients of prototypes

		d <- data.frame(
			me = centers, 
			time.point = max.col(t(x$landscape$TMC[, centers, drop = FALSE])),	# the time points of each active prototypes
			parent = 0,
			child = 0 
		)
		d[d[, 'time.point'] == 1, 'parent'] <- NA	# the cells at the first time point has no parent
		d[d[, 'time.point'] == nrow(x$landscape$TMC), 'child'] <- NA	# the cells at the last time point has no children
		for (i in 1:nrow(d)){
			if (!is.na(d[i, 'parent']) && d[i, 'parent'] == 0){
				h <- which(x$landscape$TMC[d[i, 'time.point'] - 1, ] & is.center)	# active prototypes at the previous layer
				parent <- h[knnx.index(t(Theta[, h]), t(Theta[, d[i, 'me'], drop = FALSE]), k = 1)[, 1]]	# the nearest convolving prototypes at current layer
				d[i, 'parent'] <- parent
			}
			if (!is.na(d[i, 'child']) && d[i, 'child'] == 0){
				h <- which(x$landscape$TMC[d[i, 'time.point'] + 1, ] & is.center)	# active prototypes at the next layer
				child <- h[knnx.index(t(Theta[, h]), t(Theta[, d[i, 'me'], drop = FALSE]), k = 1)[, 1]]	# the nearest convolving prototypes at current layer
				d[i, 'child'] <- child
			}
		}

		g <- rbind(
			data.frame(from = d[, 'me'], to = d[, 'child']),
			data.frame(from = d[, 'parent'], to = d[, 'me'])
		)
		g <- unique(g[!is.na(g[, 'from']) & !is.na(g[, 'to']), ])
		g <- do.call('rbind', lapply(1:nrow(g), function(i){
			vs <- as.vector(get.shortest.paths(x$landscape$graph, from = g[i, 'from'], to = g[i, 'to'])$vpath[[1]])
			cbind(from = vs[1:(length(vs) - 1)], to = vs[2:length(vs)])
		}))
		G <- sparseMatrix(i = g[, 'from'], j = g[, 'to'], x = sqrt(colSums((Theta[, g[, 'from']] - Theta[, g[, 'to']])^2)), dims = c(x$landscape$H.prototype, x$landscape$H.prototype))
		G <- (G + t(G)) / 2
		G <- graph.adjacency(G, mode = 'undirected', weighted = TRUE)
		roots <- d[is.na(d[, 'parent']), 'me']
		leaves <- d[is.na(d[, 'child']), 'me']
		roots <- roots[apply(distances(G, roots, leaves), 2, which.min)]
		x$landscape$paths <- lapply(1:length(leaves), function(i) as.vector(get.shortest.paths(G, from = roots[i], to = leaves[i])$vpath[[1]]))
		x$branching <- do.call('cbind', lapply(x$landscape$paths, function(i) mem %in% i))
		is.orphan <- rowSums(x$branching > 0) == 0
		if (any(is.orphan))
			x$branching[is.orphan, ] <- NA
	}
	x

} # end of trajectory


#' fastmds
#' 
#' A sampling-based fastMDS algorithm based on Yang et al.
#'
#' @param X The input gene by cell matrix
#' @param K Number of dimensions
#' @param f a vector indicating the pre-defined groups of input data
#' @param scale whether or not rescale the initialized metagene coefficients (default: FALSE)
#' @param control control parameters
#'
fastmds <- function(X, K, f = NULL, scale = FALSE, control = NULL){

	M <- ncol(X)
	if (any(table(f) <= K))
		stop('pre-defined group size is smaller than K')

	control <- set.control(control)

	s <- split.data(M, f = f, control = control)	# split the data

	if (M < control$max.size.dist.matrix && is.null(f)){
		V <- tryCatch({
			t(cmdscale(as.dist(rdist(t(X))), eig = FALSE, k = K))
		}, error = function(e){
			stop('cmdscale failed')
		})
	}else{
		control$max.size.per.batch <- min(control$max.size.dist.matrix, control$max.size.per.batch)
		V.list <- mclapply(1:s$n.batch, function(b){
			t(cmdscale(as.dist(rdist(t(X[, s$groups == b]))), eig = FALSE, k = K))
		}, mc.cores = control$mc.cores)	# MDS of cells within each batch
		m <- lapply(1:s$n.batch, function(b) sample(1:s$size[b], max(2, min(s$size[b], ceiling(control$max.size.dist.matrix / s$n.batch)))))	# for sampling a subset of cells from each batch
		m.align <- lapply(1:s$n.batch, function(b) which(s$groups == b)[m[[b]]])	# convert the local index to global index
		V.align <- t(cmdscale(as.dist(rdist(t(X[, unlist(m.align)]))), eig = FALSE, k = K))	# compute the MDS of sampled cells
		V.align <- lapply(split(1:ncol(V.align), list(rep(1:s$n.batch, sapply(m, length)))), function(i) V.align[, i, drop = FALSE])
		V.list <- mclapply(1:s$n.batch, function(b){	# map the local MDS to global MDS
			s <- svd(V.list[[b]][, m[[b]]])
			V.align[[b]] %*% s$v %*% diag(1 / s$d) %*% t(s$u) %*% V.list[[b]]
		}, mc.cores = control$mc.cores)
		V <- matrix(0, K, M)
		for (b in 1:s[['n.batch']]){
			V[, s[['groups']] == b] <- V.list[[b]]
		}
	}
	if (scale) 
		V <- scale(V)
	V <- matrix(c(V), nrow(V), ncol(V), dimnames = list(NULL, colnames(X)))
	V
} # end of fastmds



#' gtm
#'
#' Update the prototype landscape
#'
gtm <- function(V = NULL, landscape = NULL, CT = NULL, method = NULL, update.beta = TRUE, control = NULL){

	M <- ncol(V)	# number of cells

	if (landscape$type == 'plate')
		method <- 'all'
	
	if (is.null(method))
		method <- 'all'

	if (is.null(update.beta))
		update.beta <- TRUE

	if (!update.beta)
		V <- scale(V)

	# decide the order of the cells for updating the landscape
	if (method == 'forward'){
		P <- do.call('rbind', lapply(1:ncol(CT), function(t) Matrix::rowSums(CT[, 1:t, drop = FALSE]) > 0))
	}else if (method %in% c('backward')){
		P <- do.call('rbind', lapply(ncol(CT):1, function(t) Matrix::rowSums(CT[, t:ncol(CT), drop = FALSE]) > 0))
	}else if (method == 'all'){
		P <- matrix(TRUE, nrow = 1, ncol = M)
	}else
		stop(sprintf('unknown method: %s', method))

	for (i in 1:nrow(P)){
		if (i == 1)
			mem <- NULL
		else{
			memp <- max.col(mf$Z)
			memp[apply(mf$Z, 1, max) < landscape$assignment.threshold] <- NA
			mem <- rep(NA, sum(P[i, ]))
			mem[P[i - 1, ]] <- memp
		}
		mf <- gtm.core(V = V[, P[i, ]], landscape = landscape, CT = CT[P[i, ], ], mem = mem[P[i, ]], update.beta = update.beta, control = control)
		landscape <- mf$landscape
	}
	mf

} # end of update.landscape


#' Get a time table that consider the pseodotime information
pseudotime.table <- function(X, n = NA, time.table = NULL, K = 3, direction = 'center', control = NULL){

	eps <- .Machine$double.eps
	M <- ncol(X)

	if (is.na(n) || n <= 1)
		stop('n must be greater than 1')

	V <- fastmds(scale(log(X + 1)), K = K, control = control)
	V <- t(apply(V, 1, function(v) (v - min(v)) / (max(v) - min(v))))
	pt.order <- colMeans(abs(V - 0.5))
	pt.order <- (pt.order - min(pt.order)) / (max(pt.order) - min(pt.order))

	time.order <- max.col(time.table)
	time.order <- time.order / ncol(time.table)
	cc <- cor(pt.order, time.order, method = 'kendall')
	if (cc < 0)
		pt.order <- 1 - pt.order
	cc <- abs(cc)
	w <- min((cc / 0.2)^2, 1)	# weight for time order
	
	od <- (1 - w) * pt.order + w * time.order
	od <- (od - min(od)) / (max(od) - min(od))

	breaks <- quantile(od, seq(0, 1, length.out = n + 1))
	breaks[1] <- breaks[1] - 1
	breaks[length(breaks)] <- breaks[length(breaks)] + 1
	sparseMatrix(i = 1:M, j = as.numeric(cut(od, breaks)), dims = c(M, n))

} # end of pseudo.time.table



plot.landscape <- function(x, ...){

	if (x$type == 'plate'){

		H.prototype <- x$H.prototype	# number of prototypes
		coord <- x$Y.prototype
		Yp <- x$Yp.prototype
		n.prototype <- x$n.prototype

		lim <- c(-max(Yp[, 'r']), max(Yp[, 'r']))	# limit of the 2D space
		plot(NA, xlim = lim * 1.1, ylim = lim * 1.1, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', asp = 1)

		rs <- unique(Yp[, 'r'])
		col.time <- colorpanel(length(rs), low = 'gray100', mid = 'gray70', high = 'gray40')
		for (i in length(rs):1){
			draw.circle(0, 0, rs[i], nv = 100, border = NA, col = col.time[i])
		}

	}
} # end of plot.landscape

