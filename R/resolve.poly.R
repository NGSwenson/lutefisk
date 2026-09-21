#' Resolve polytomies with informative branch lengths
#'
#' \code{\link[ape]{multi2di}} resolves a polytomy into a random sequence of
#' bifurcations but places every new node at zero distance from the original
#' polytomous node, so the resolved tree has the same phylogenetic
#' variance-covariance matrix and the same pairwise distances as the
#' unresolved one. This function resolves the topology randomly in the same
#' way and then gives the new nodes positions along the branches that lead to
#' the original descendants, so that the resolutions differ in branch length
#' as well as in topology. Node positions are drawn in one of three ways and
#' the tree remains ultrametric if it was ultrametric to begin with.
#'
#' Each new node sits between its parent (the polytomous node or another new
#' node above it) and the shallowest of the original nodes below it, and its
#' position is expressed as a fraction \eqn{u} of that interval measured from
#' the parent. With \code{method = "equidistant"} the new nodes are spaced
#' evenly along the longest chain of new nodes below the polytomy, which is
#' the deterministic choice and the one to use when a single resolved tree
#' is wanted. With \code{method = "random"} the fraction is drawn from
#' \code{positions}: by default a uniform draw on (0, 1), or a vector of
#' fractions to sample from, with optional sampling weights \code{prob}, or
#' a function of one argument (the number of draws) returning fractions.
#' With \code{method = "directed"} the fraction is drawn from a beta
#' distribution with shape parameters \code{shape}, which lets the user
#' push the new nodes toward the polytomous node (for example
#' \code{shape = c(1, 4)}, appropriate when the polytomy is thought to be a
#' rapid radiation) or toward the descendants (\code{shape = c(4, 1)}).
#'
#' @param my.phylo A \code{phylo} object with branch lengths.
#' @param method \code{"equidistant"}, \code{"random"} or \code{"directed"}.
#' @param n Number of resolved trees to return.
#' @param positions For \code{method = "random"}: \code{NULL} (uniform), a
#'   numeric vector of fractions in (0, 1) to sample from, or a function.
#' @param prob Optional sampling weights when \code{positions} is a vector.
#' @param shape Two positive shape parameters for \code{method = "directed"}.
#' @param random.topology Resolve the topology at random (the default) or in
#'   the fixed order that \code{\link[ape]{multi2di}} uses when
#'   \code{random = FALSE}.
#' @return A \code{phylo} object when \code{n = 1}, otherwise a
#'   \code{multiPhylo} object of length \code{n}.
#' @seealso \code{\link{perturb.branches}} for uncertainty in the branch
#'   lengths of an already bifurcating tree.
#' @examples
#' data(example.poly)
#' r1 <- resolve.poly(example.poly, method = "equidistant")
#' r2 <- resolve.poly(example.poly, method = "random", n = 100)
#' r3 <- resolve.poly(example.poly, method = "directed", shape = c(1, 4), n = 100)
#' ape::is.ultrametric(r1)
#' range(sapply(r2, function(x) cophenetic(x)["sp1", "sp2"]))
#' @export
resolve.poly <- function(my.phylo, method = c("equidistant", "random", "directed"),
	n = 1, positions = NULL, prob = NULL, shape = c(1, 1), random.topology = TRUE) {
	method <- match.arg(method)
	if (is.null(my.phylo$edge.length)) stop("my.phylo must have branch lengths")
	if (ape::is.binary(my.phylo)) {
		warning("my.phylo has no polytomies; returning it unchanged")
		if (n == 1) return(my.phylo)
		out <- rep(list(my.phylo), n); class(out) <- "multiPhylo"; return(out)
	}
	draw <- switch(method,
		equidistant = NULL,
		random = {
			if (is.null(positions)) {
				function(k) stats::runif(k)
			} else if (is.function(positions)) {
				positions
			} else {
				if (any(positions <= 0 | positions >= 1)) {
					stop("positions must lie strictly between 0 and 1")
				}
				function(k) sample(positions, k, replace = TRUE, prob = prob)
			}
		},
		directed = function(k) stats::rbeta(k, shape[1], shape[2]))

	## Tag the original nodes so that the new ones can be told apart after
	## multi2di() renumbers everything.
	tagged <- my.phylo
	tagged$node.label <- paste0(".orig", seq_len(my.phylo$Nnode))
	had.labels <- !is.null(my.phylo$node.label)

	one <- function() {
		di <- ape::multi2di(tagged, random = random.topology)
		n.tip <- length(di$tip.label)
		is.new <- c(rep(FALSE, n.tip), di$node.label == "")
		## Heights of every node measured from the root.
		h <- ape::node.depth.edgelength(di)
		parent <- integer(n.tip + di$Nnode)
		parent[di$edge[, 2]] <- di$edge[, 1]
		kids <- split(di$edge[, 2], di$edge[, 1])
		## Upper bound for each new node: the shallowest original node below
		## it. Computed bottom-up so that chains of new nodes resolve.
		upper <- h
		chain <- integer(length(h))
		new.nodes <- which(is.new)
		## process deepest first (postorder over new nodes by height is not
		## reliable since all have equal height; use descendant structure)
		remaining <- new.nodes
		done <- logical(length(h))
		done[!is.new] <- TRUE
		while (length(remaining) > 0) {
			ready <- remaining[sapply(remaining, function(v) all(done[kids[[as.character(v)]]]))]
			for (v in ready) {
				k <- kids[[as.character(v)]]
				upper[v] <- min(upper[k])
				chain[v] <- 1 + max(chain[k])
				done[v] <- TRUE
			}
			remaining <- setdiff(remaining, ready)
		}
		## Assign heights top-down: the parent of a new node is either the
		## original polytomous node or a new node already assigned.
		ord <- new.nodes[order(-chain[new.nodes])]
		for (v in ord) {
			lower <- h[parent[v]]
			u <- if (is.null(draw)) 1 / (chain[v] + 1) else draw(1)
			h[v] <- lower + u * (upper[v] - lower)
		}
		di$edge.length <- h[di$edge[, 2]] - h[di$edge[, 1]]
		if (had.labels) {
			lab <- rep("", di$Nnode)
			orig <- di$node.label != ""
			lab[orig] <- my.phylo$node.label[as.integer(sub(".orig", "", di$node.label[orig]))]
			di$node.label <- lab
		} else {
			di$node.label <- NULL
		}
		di
	}
	if (n == 1) return(one())
	out <- replicate(n, one(), simplify = FALSE)
	class(out) <- "multiPhylo"
	out
}

#' Perturb branch lengths to represent dating uncertainty
#'
#' Multiplies every branch length by an independent lognormal factor and,
#' optionally, rescales the result so that the tree keeps its original
#' depth. This is a crude substitute for a posterior sample of trees when
#' only a single dated tree is available: it generates trees whose
#' phylogenetic variance-covariance matrices differ, which is what
#' VCV-based methods such as phylogenetic signal and PGLS need in order to
#' show any sensitivity at all to the tree.
#'
#' @param my.phylo A \code{phylo} object with branch lengths.
#' @param sdlog Standard deviation of the log of the multiplicative factor.
#' @param n Number of trees to return.
#' @param keep.depth Rescale each tree so its maximum root-to-tip distance
#'   equals that of the original tree.
#' @return A \code{phylo} object when \code{n = 1}, otherwise a
#'   \code{multiPhylo} object.
#' @examples
#' data(my.phylo)
#' trees <- perturb.branches(my.phylo, sdlog = 0.3, n = 10)
#' sapply(trees, function(x) cophenetic(x)["a", "b"])
#' @export
perturb.branches <- function(my.phylo, sdlog = 0.3, n = 1, keep.depth = TRUE) {
	depth <- max(ape::node.depth.edgelength(my.phylo))
	one <- function() {
		tr <- my.phylo
		tr$edge.length <- tr$edge.length * stats::rlnorm(length(tr$edge.length), 0, sdlog)
		if (keep.depth) {
			tr$edge.length <- tr$edge.length * depth / max(ape::node.depth.edgelength(tr))
		}
		tr
	}
	if (n == 1) return(one())
	out <- replicate(n, one(), simplify = FALSE)
	class(out) <- "multiPhylo"
	out
}
