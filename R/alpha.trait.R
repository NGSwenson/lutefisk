#' Moments of trait distributions within communities
#'
#' The mean, standard deviation, skewness and kurtosis of each trait across
#' the species present in each community (\code{com.trait.moments}), or the
#' abundance-weighted mean and standard deviation
#' (\code{com.trait.weighted}). The first edition looped over traits and
#' applied over communities; here every moment for every trait and every
#' community comes from products of the community matrix with powers of the
#' trait matrix, so the functions scale to thousands of communities.
#'
#' Skewness is the third central moment divided by the cube of the
#' (population) standard deviation and kurtosis is the fourth central
#' moment divided by the squared variance minus three, as in the first
#' edition. Weighted standard deviations use \code{\link{weighted.sd}}.
#'
#' @param my.sample A community data matrix with communities as rows and
#'   species as columns.
#' @param traits A numeric matrix or data frame of traits with species as
#'   row names. Species with missing values for a trait are dropped from
#'   the calculation of that trait.
#' @return A matrix with one row per community and one column per moment
#'   and trait combination.
#' @examples
#' data(lefse.sample, lefse.traits)
#' com.trait.moments(lefse.sample, lefse.traits)
#' com.trait.weighted(lefse.sample, lefse.traits)
#' @export
com.trait.moments <- function(my.sample, traits) {
	pa <- to.pa(my.sample)
	tr <- as.matrix(traits)[colnames(pa), , drop = FALSE]
	storage.mode(tr) <- "numeric"
	## a separate presence matrix per trait so that NA trait values are
	## simply not counted
	m1 <- m2 <- m3 <- m4 <- matrix(NA_real_, nrow(pa), ncol(tr))
	for (j in seq_len(ncol(tr))) {
		ok <- !is.na(tr[, j])
		p <- pa[, ok, drop = FALSE]
		t1 <- tr[ok, j]
		n <- rowSums(p)
		mu <- as.vector(p %*% t1) / n
		## central moments: for community i, sum over present species of
		## (t - mu_i)^k. Expand the power so that one product suffices.
		s1 <- as.vector(p %*% t1)
		s2 <- as.vector(p %*% t1^2)
		s3 <- as.vector(p %*% t1^3)
		s4 <- as.vector(p %*% t1^4)
		c2 <- (s2 - n * mu^2) / n
		c3 <- (s3 - 3 * mu * s2 + 3 * mu^2 * s1 - n * mu^3) / n
		c4 <- (s4 - 4 * mu * s3 + 6 * mu^2 * s2 - 4 * mu^3 * s1 + n * mu^4) / n
		m1[, j] <- mu
		m2[, j] <- sqrt(c2 * n / (n - 1))
		m3[, j] <- c3 / (c2 * n / (n - 1))^(3 / 2)
		m4[, j] <- c4 / (c2 * n / (n - 1))^2 - 3
	}
	out <- cbind(m1, m2, m3, m4)
	colnames(out) <- paste(rep(c("mean", "sd", "skew", "kurtosis"),
		each = ncol(tr)), rep(colnames(tr), 4), sep = ".")
	rownames(out) <- rownames(pa)
	out
}

#' @rdname com.trait.moments
#' @export
com.trait.weighted <- function(my.sample, traits) {
	ra <- to.ra(my.sample)
	tr <- as.matrix(traits)[colnames(ra), , drop = FALSE]
	storage.mode(tr) <- "numeric"
	wm <- wsd <- matrix(NA_real_, nrow(ra), ncol(tr))
	for (j in seq_len(ncol(tr))) {
		ok <- !is.na(tr[, j])
		w <- ra[, ok, drop = FALSE]
		w <- w / rowSums(w)
		t1 <- tr[ok, j]
		mu <- as.vector(w %*% t1)
		wm[, j] <- mu
		## reliability-weights variance, as in SDMTools::wt.var
		ss <- as.vector(w %*% t1^2) - mu^2
		wsd[, j] <- sqrt(ss * (1 / (1 - rowSums(w^2))))
	}
	out <- cbind(wm, wsd)
	colnames(out) <- paste(rep(c("mean", "sd"), each = ncol(tr)),
		rep(colnames(tr), 2), sep = ".")
	rownames(out) <- rownames(ra)
	out
}

#' Trait ranges within communities
#'
#' The range (maximum minus minimum) of each trait across the species
#' present in each community.
#'
#' @inheritParams com.trait.moments
#' @return A matrix with one row per community and one column per trait.
#' @examples
#' data(lefse.sample, lefse.traits)
#' trait.range(lefse.sample, lefse.traits)
#' @export
trait.range <- function(my.sample, traits) {
	pa <- to.pa(my.sample)
	tr <- as.matrix(traits)[colnames(pa), , drop = FALSE]
	storage.mode(tr) <- "numeric"
	out <- t(apply(pa, 1, function(p) {
		sub <- tr[p > 0, , drop = FALSE]
		apply(sub, 2, max, na.rm = TRUE) - apply(sub, 2, min, na.rm = TRUE)
	}))
	dimnames(out) <- list(rownames(pa), colnames(tr))
	out
}

#' Convex hull volume of communities in trait space
#'
#' The volume of the convex hull enclosing the species present in each
#' community, in the trait space supplied. This is the uncorrected
#' functional richness of Cornwell et al. (2006). The \code{dbFD} function in
#' \pkg{FD} computes its FRic on principal coordinate axes derived from a
#' Gower distance and rescales it; this function works directly on the
#' trait matrix given, so the two agree only when the same space is used.
#'
#' @inheritParams com.trait.moments
#' @return A named numeric vector of hull volumes. Communities with fewer
#'   present species than \code{ncol(traits) + 1} return \code{NA}.
#' @examples
#' data(lefse.sample, lefse.traits)
#' hull.volume(lefse.sample, scale(lefse.traits))
#' @export
hull.volume <- function(my.sample, traits) {
	if (!requireNamespace("geometry", quietly = TRUE)) {
		stop("the geometry package is required for hull.volume")
	}
	pa <- to.pa(my.sample)
	tr <- as.matrix(traits)[colnames(pa), , drop = FALSE]
	storage.mode(tr) <- "numeric"
	out <- apply(pa, 1, function(p) {
		sub <- tr[p > 0, , drop = FALSE]
		sub <- sub[stats::complete.cases(sub), , drop = FALSE]
		if (nrow(sub) <= ncol(sub)) return(NA_real_)
		geometry::convhulln(sub, options = "FA")$vol
	})
	names(out) <- rownames(pa)
	out
}
