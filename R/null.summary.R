#' Summarize a null model
#'
#' Standardized effect sizes, robust standardized effect sizes (using the
#' median and median absolute deviation of the null distribution), ranks and
#' one-tailed probabilities for observed values against a matrix of null
#' values.
#'
#' @param obs A numeric vector with one observed value per community.
#' @param nulls A numeric matrix with one row per community and one column
#'   per randomization.
#' @return A data frame with one row per community.
#' @examples
#' data(lefse.sample, lefse.phylo)
#' d <- cophenetic(lefse.phylo)
#' obs <- Fmpd(d, lefse.sample)
#' set.seed(1)
#' nulls <- replicate(199, Fmpd(shuffle.dist(d), lefse.sample))
#' null.summary(obs, nulls)
#' @export
null.summary <- function(obs, nulls) {
	nulls <- as.matrix(nulls)
	n.rand <- ncol(nulls)
	null.mean <- rowMeans(nulls, na.rm = TRUE)
	null.sd <- apply(nulls, 1, stats::sd, na.rm = TRUE)
	null.med <- apply(nulls, 1, stats::median, na.rm = TRUE)
	null.mad <- apply(nulls, 1, stats::mad, na.rm = TRUE)
	obs.rank <- apply(cbind(obs, nulls), 1, function(x) rank(x)[1])
	data.frame(obs = obs, null.mean = null.mean, null.sd = null.sd,
		ses = (obs - null.mean) / null.sd,
		ses.robust = (obs - null.med) / null.mad,
		rank = obs.rank, p.lower = obs.rank / (n.rand + 1),
		p.upper = 1 - (obs.rank - 1) / (n.rand + 1), n.rand = n.rand)
}
