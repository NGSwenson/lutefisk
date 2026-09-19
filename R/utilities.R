#' Tip-by-branch membership matrix for a phylogeny
#'
#' Build a binary matrix with one row per tip and one column per branch
#' (edge) of a phylogeny, in which an entry is 1 when the tip descends from
#' that branch. Most tree-based diversity metrics in this package are
#' computed as products of a community data matrix with this matrix, which
#' replaces the per-branch loops used in the first edition of the book.
#'
#' @param my.phylo A \code{phylo} object.
#' @return An integer matrix with \code{Ntip(my.phylo)} rows named by tip
#'   label and \code{nrow(my.phylo$edge)} columns, one per row of
#'   \code{my.phylo$edge}.
#' @examples
#' data(PD.example.phylo)
#' m <- branch.membership(PD.example.phylo)
#' dim(m)
#' ## the number of tips subtended by each branch
#' colSums(m)
#' @export
branch.membership <- function(my.phylo) {
	n.tip <- length(my.phylo$tip.label)
	edge <- my.phylo$edge
	## parent[v] is the parent node of node v (0 for the root).
	parent <- integer(n.tip + my.phylo$Nnode)
	parent[edge[, 2]] <- edge[, 1]
	## edge.of[v] is the row of the edge matrix ending in node v.
	edge.of <- integer(n.tip + my.phylo$Nnode)
	edge.of[edge[, 2]] <- seq_len(nrow(edge))
	out <- matrix(0L, nrow = n.tip, ncol = nrow(edge),
		dimnames = list(my.phylo$tip.label, NULL))
	## Walk from each tip to the root and mark every branch on the path.
	for (i in seq_len(n.tip)) {
		v <- i
		while (parent[v] != 0) {
			out[i, edge.of[v]] <- 1L
			v <- parent[v]
		}
	}
	out
}

#' Tip-by-node membership matrix for a phylogeny
#'
#' As \code{\link{branch.membership}} but with one column per internal node
#' rather than per branch. Used by \code{\link{lineage.hyper}}.
#'
#' @inheritParams branch.membership
#' @return An integer matrix with one row per tip and one column per
#'   internal node, columns named by node number.
#' @export
node.membership <- function(my.phylo) {
	n.tip <- length(my.phylo$tip.label)
	m <- branch.membership(my.phylo)
	## Column j of the branch matrix ends in node edge[j, 2]; keep those that
	## end in internal nodes, plus a column of ones for the root.
	internal <- my.phylo$edge[, 2] > n.tip
	out <- cbind(1L, m[, internal, drop = FALSE])
	colnames(out) <- c(n.tip + 1, my.phylo$edge[internal, 2])
	out[, order(as.integer(colnames(out))), drop = FALSE]
}

#' Convert a community data matrix to presence-absence or relative abundance
#'
#' @param my.sample A community data matrix with communities as rows and
#'   species as columns.
#' @return A numeric matrix of the same dimensions.
#' @examples
#' data(lefse.sample)
#' to.pa(lefse.sample)
#' rowSums(to.ra(lefse.sample))
#' @export
to.pa <- function(my.sample) {
	x <- as.matrix(my.sample)
	storage.mode(x) <- "numeric"
	(x > 0) * 1
}

#' @rdname to.pa
#' @export
to.ra <- function(my.sample) {
	x <- as.matrix(my.sample)
	storage.mode(x) <- "numeric"
	x / rowSums(x)
}

#' Weighted standard deviation
#'
#' A weighted standard deviation with the small-sample correction for
#' reliability weights, matching the \code{wt.sd} function of the retired
#' \pkg{SDMTools} package that the first edition of the book relied on.
#'
#' @param x A numeric vector.
#' @param w A vector of non-negative weights of the same length.
#' @return The weighted standard deviation.
#' @examples
#' weighted.sd(1:5, c(1, 1, 1, 1, 10))
#' @export
weighted.sd <- function(x, w) {
	ok <- is.finite(x * w)
	x <- x[ok]
	w <- w[ok]
	xbar <- sum(w * x) / sum(w)
	v <- sum(w * (x - xbar)^2) * (sum(w) / (sum(w)^2 - sum(w^2)))
	sqrt(v)
}

## Check that a distance matrix covers every species in a sample and return
## it as a numeric matrix in sample column order.
match.dist <- function(dist.mat, my.sample) {
	d <- as.matrix(dist.mat)
	spp <- colnames(my.sample)
	missing <- setdiff(spp, rownames(d))
	if (length(missing) > 0) {
		stop("species missing from the distance matrix: ",
			paste(missing, collapse = ", "))
	}
	d[spp, spp]
}
