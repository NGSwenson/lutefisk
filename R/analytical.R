#' Exact null moments of Faith's index and mean pairwise distance
#'
#' Under the null model that a community of \eqn{n} species is a uniform
#' random draw of \eqn{n} tips from a pool of \eqn{N}, the expectation and
#' variance of Faith's phylogenetic diversity and of the mean pairwise
#' distance have closed forms, so a standardized effect size can be computed
#' exactly and instantly rather than estimated by randomization. These are
#' the quantities that the \code{PhyloMeasures} package
#' (Tsirogiannis and Sandel 2016) computed, and the derivations follow
#' Tsirogiannis et al. (2012). They apply to the name-shuffling
#' (\code{"taxa.labels"}) null model without abundance weighting; for
#' constrained null models, and for abundance-weighted metrics, there is no
#' closed form and one must randomize.
#'
#' For Faith's index, each branch \eqn{b} of length \eqn{L_b} subtending
#' \eqn{t_b} tips is included in the community's subtree with probability
#' \eqn{p_b = 1 - {N - t_b \choose n} / {N \choose n}}, so
#' \eqn{E[PD] = \sum_b L_b p_b}. The variance requires the probability that
#' two branches are both included, which depends on the number of tips
#' under their union, and the whole calculation is a few matrix products on
#' the branch membership matrix. For the mean pairwise distance the
#' expectation is the mean of all pairwise distances in the pool, whatever
#' \eqn{n}, and the variance is a function of the sums of the distances, of
#' their squares, and of the squared row sums of the distance matrix.
#'
#' @param my.phylo A \code{phylo} object for the species pool.
#' @param n A vector of community sizes (numbers of species).
#' @param dist.mat For \code{expected.mpd}, the pairwise distance matrix of
#'   the pool (any symmetric distance matrix, phylogenetic or functional).
#' @return A data frame with one row per value of \code{n} giving
#'   \code{n}, the expected value and the standard deviation.
#' @references Tsirogiannis, C., Sandel, B. and Cheliotis, D. (2012)
#'   Efficient computation of popular phylogenetic tree measures. Lecture
#'   Notes in Computer Science 7534, 30-43. Tsirogiannis, C. and Sandel, B.
#'   (2016) PhyloMeasures: a package for computing phylogenetic biodiversity
#'   measures and their statistical moments. Ecography 39, 709-714.
#' @examples
#' data(null.pd.example.phylo)
#' expected.pd(null.pd.example.phylo, n = c(5, 10, 20))
#' expected.mpd(cophenetic(null.pd.example.phylo), n = c(5, 10, 20))
#' @export
expected.pd <- function(my.phylo, n) {
	m <- branch.membership(my.phylo)
	N <- nrow(m)
	L <- my.phylo$edge.length
	t.b <- colSums(m)
	## number of tips under the union of each pair of branches
	shared <- crossprod(m)
	union <- outer(t.b, t.b, "+") - shared
	out <- t(sapply(n, function(k) {
		if (k > N) stop("n cannot exceed the number of tips")
		## probability that a set of s tips is missed entirely by a draw of k
		miss <- function(s) exp(lchoose(N - s, k) - lchoose(N, k))
		p.b <- 1 - miss(t.b)
		e <- sum(L * p.b)
		p.both <- 1 - miss(t.b) - rep(miss(t.b), each = length(t.b)) + miss(union)
		dim(p.both) <- dim(union)
		cov.b <- p.both - outer(p.b, p.b)
		v <- as.numeric(t(L) %*% cov.b %*% L)
		c(n = k, expected = e, sd = sqrt(max(v, 0)))
	}))
	as.data.frame(out)
}

#' @rdname expected.pd
#' @export
expected.mpd <- function(dist.mat, n) {
	d <- as.matrix(dist.mat)
	N <- nrow(d)
	diag(d) <- 0
	T1 <- sum(d[upper.tri(d)])
	Q <- sum(d[upper.tri(d)]^2)
	R2 <- sum(rowSums(d)^2)
	out <- t(sapply(n, function(k) {
		if (k > N) stop("n cannot exceed the number of species")
		if (k < 2) return(c(n = k, expected = NA, sd = NA))
		p2 <- k * (k - 1) / (N * (N - 1))
		p3 <- k * (k - 1) * (k - 2) / (N * (N - 1) * (N - 2))
		p4 <- if (N > 3) k * (k - 1) * (k - 2) * (k - 3) / (N * (N - 1) * (N - 2) * (N - 3)) else 0
		## X = sum of distances over pairs in the sample
		ex <- p2 * T1
		ex2 <- p2 * Q + p3 * (R2 - 2 * Q) + p4 * (T1^2 + Q - R2)
		pairs <- k * (k - 1) / 2
		c(n = k, expected = ex / pairs, sd = sqrt(max(ex2 - ex^2, 0)) / pairs)
	}))
	as.data.frame(out)
}

#' Analytical standardized effect sizes
#'
#' Computes the observed Faith's index or mean pairwise distance for every
#' community and standardizes it by the exact null mean and standard
#' deviation from \code{\link{expected.pd}} or \code{\link{expected.mpd}},
#' together with a two-tailed probability from the normal approximation.
#'
#' @inheritParams expected.pd
#' @param my.sample A community data matrix; only presence is used.
#' @param metric \code{"pd"} or \code{"mpd"}.
#' @return A data frame with one row per community.
#' @examples
#' data(null.pd.example.phylo)
#' set.seed(1)
#' cdm <- t(replicate(4, as.integer(runif(40) < 0.3)))
#' colnames(cdm) <- null.pd.example.phylo$tip.label
#' ses.analytical(null.pd.example.phylo, cdm, metric = "pd")
#' @export
ses.analytical <- function(my.phylo, my.sample, metric = c("pd", "mpd"),
	dist.mat = NULL) {
	metric <- match.arg(metric)
	pa <- to.pa(my.sample)
	pa <- pa[, my.phylo$tip.label, drop = FALSE]
	n <- rowSums(pa)
	sizes <- sort(unique(n))
	if (metric == "pd") {
		obs <- faith.pd(my.phylo, pa)
		mom <- expected.pd(my.phylo, sizes)
	} else {
		if (is.null(dist.mat)) dist.mat <- stats::cophenetic(my.phylo)
		obs <- Fmpd(dist.mat, pa)
		mom <- expected.mpd(dist.mat, sizes)
	}
	idx <- match(n, mom$n)
	z <- (obs - mom$expected[idx]) / mom$sd[idx]
	data.frame(n = n, obs = obs, expected = mom$expected[idx],
		sd = mom$sd[idx], ses = z, p = 2 * stats::pnorm(-abs(z)),
		row.names = rownames(pa))
}
