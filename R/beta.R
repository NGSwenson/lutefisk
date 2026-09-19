#' Pairwise phylogenetic or functional beta diversity between communities
#'
#' Mean pairwise distance (\code{new.comdist}, \code{new.comdist.prime}) and
#' mean nearest neighbor distance (\code{new.comdistnn},
#' \code{new.comdistnn.prime}) between the species of every pair of
#' communities, unweighted or weighted by relative abundance (the
#' \code{.prime} versions). The first edition computed each pair inside a
#' nested \code{lapply}; here the whole community-by-community matrix comes
#' from one or two matrix products.
#'
#' For presence-absence matrix \eqn{P} and distance matrix \eqn{D}, the sum
#' of distances between all species of community \eqn{i} and all species of
#' community \eqn{j} is element \eqn{(i, j)} of \eqn{P D P'}, so
#' \code{new.comdist} is that matrix divided by the outer product of the
#' species richnesses. The abundance-weighted version substitutes relative
#' abundances for \eqn{P}, which removes the need for a divisor.
#'
#' For the nearest neighbor versions, the distance from each species
#' \eqn{s} to its nearest neighbor in community \eqn{j} is computed once for
#' every species and community, giving a species-by-community matrix
#' \eqn{N}. Then \eqn{P N} sums, for each pair \eqn{(i, j)}, the nearest
#' neighbor distances of the species of \eqn{i} into \eqn{j}, and the
#' symmetric mean of that matrix and its transpose divided by the summed
#' richness is \code{new.comdistnn}. The weighted version uses relative
#' abundances and divides by two.
#'
#' @param my.sample A community data matrix with communities as rows and
#'   species as columns.
#' @param my.dist.mat A species-by-species distance matrix.
#' @return A community-by-community matrix. For the pairwise versions the
#'   diagonal is the within-community mean pairwise distance including
#'   self comparisons; for the nearest neighbor versions it is zero.
#' @examples
#' data(beta.example.sample, beta.example.phylo)
#' d <- cophenetic(beta.example.phylo)
#' new.comdist(beta.example.sample, d)
#' new.comdistnn.prime(beta.example.sample, d)
#' @export
new.comdist <- function(my.sample, my.dist.mat) {
	pa <- to.pa(my.sample)
	d <- match.dist(my.dist.mat, pa)
	n <- rowSums(pa)
	(pa %*% d %*% t(pa)) / outer(n, n)
}

#' @rdname new.comdist
#' @export
new.comdist.prime <- function(my.sample, my.dist.mat) {
	ra <- to.ra(my.sample)
	d <- match.dist(my.dist.mat, ra)
	ra %*% d %*% t(ra)
}

## Species-by-community matrix of nearest neighbor distances: entry [s, j]
## is the distance from species s to the closest species present in
## community j. A species shared by both communities is its own nearest
## neighbor at distance zero, as in picante::comdistnt, so the diagonal of
## the resulting beta diversity matrix is zero.
nn.between <- function(dist.mat, pa) {
	d <- match.dist(dist.mat, pa)
	out <- apply(pa, 1, function(p) {
		present <- p > 0
		apply(d[, present, drop = FALSE], 1, min)
	})
	dimnames(out) <- list(colnames(pa), rownames(pa))
	out
}

#' @rdname new.comdist
#' @export
new.comdistnn <- function(my.sample, my.dist.mat) {
	pa <- to.pa(my.sample)
	nn <- nn.between(my.dist.mat, pa)
	n <- rowSums(pa)
	s <- pa %*% nn
	(s + t(s)) / outer(n, n, "+")
}

#' @rdname new.comdist
#' @export
new.comdistnn.prime <- function(my.sample, my.dist.mat) {
	ra <- to.ra(my.sample)
	nn <- nn.between(my.dist.mat, ra)
	s <- ra %*% nn
	(s + t(s)) / 2
}

#' Shuffle the species labels of a distance matrix
#'
#' A single randomization of a species-by-species distance matrix that
#' shuffles the species names on both rows and columns together. Applying a
#' diversity function to the shuffled matrix is the "taxa labels" or
#' "tip shuffle" null model, and calling this inside \code{replicate}
#' gives a null distribution. The first edition bundled one
#' randomization into each \code{.null} function; keeping the shuffle
#' separate lets any metric use it.
#'
#' @param dist.mat A species-by-species distance matrix.
#' @return The same matrix with row and column names permuted.
#' @examples
#' data(beta.example.sample, beta.example.phylo)
#' d <- cophenetic(beta.example.phylo)
#' set.seed(1)
#' nulls <- replicate(99, new.comdist(beta.example.sample, shuffle.dist(d)))
#' obs <- new.comdist(beta.example.sample, d)
#' ses <- (obs - apply(nulls, 1:2, mean)) / apply(nulls, 1:2, sd)
#' round(ses, 2)
#' @export
shuffle.dist <- function(dist.mat) {
	d <- as.matrix(dist.mat)
	new <- sample(rownames(d))
	dimnames(d) <- list(new, new)
	d
}

#' @rdname new.comdist
#' @export
new.comdist.null <- function(my.sample, my.dist.mat) {
	new.comdist(my.sample, shuffle.dist(my.dist.mat))
}

#' @rdname new.comdist
#' @export
new.comdist.prime.null <- function(my.sample, my.dist.mat) {
	new.comdist.prime(my.sample, shuffle.dist(my.dist.mat))
}

#' @rdname new.comdist
#' @export
new.comdistnn.null <- function(my.sample, my.dist.mat) {
	new.comdistnn(my.sample, shuffle.dist(my.dist.mat))
}

#' @rdname new.comdist
#' @export
new.comdistnn.prime.null <- function(my.sample, my.dist.mat) {
	new.comdistnn.prime(my.sample, shuffle.dist(my.dist.mat))
}
