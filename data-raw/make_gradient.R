## A larger simulated system for comparing metrics: 100 species on a
## phylogeny, three traits, and 40 communities arrayed along an environmental
## gradient with trait-based filtering, so that communities differ in
## phylogenetic and functional structure in ways worth comparing.
library(ape)
library(phytools)
library(geiger)
set.seed(2014)

n.spp <- 100
gradient.phylo <- pbtree(n = n.spp, scale = 1)
gradient.phylo$tip.label <- sprintf("sp%03d", seq_len(n.spp))

## Trait one evolves under Brownian motion and sets the environmental optimum;
## trait two is a weakly conserved trait (lambda about 0.3); trait three is
## drawn without regard to the phylogeny.
t1 <- as.vector(sim.char(gradient.phylo, par = 1, model = "BM")[, , 1])
lambda.tree <- rescale(gradient.phylo, model = "lambda", 0.3)
t2 <- as.vector(sim.char(lambda.tree, par = 1, model = "BM")[, , 1])
t3 <- rnorm(n.spp)
gradient.traits <- data.frame(trait1 = t1, trait2 = t2, trait3 = t3,
                              row.names = gradient.phylo$tip.label)

## 40 communities along a gradient; environment one sets the filter, the
## second is noise with spatial structure; coordinates on a transect
n.com <- 40
xy <- data.frame(x = seq(0, 100, length.out = n.com) + rnorm(n.com, 0, 2),
                 y = rnorm(n.com, 50, 15))
env1 <- seq(min(t1), max(t1), length.out = n.com)
env2 <- as.vector(scale(xy$y)) + rnorm(n.com, 0, 0.5)
gradient.env <- data.frame(xy, env1 = env1, env2 = env2,
                           row.names = sprintf("com%02d", seq_len(n.com)))

## Abundance: Poisson with a mean that falls off with the distance between a
## species' trait one value and the community's environment; the width of the
## filter varies among communities so that some are strongly filtered and
## others weakly, and a fraction of species are filtered on trait three
## instead, which is where overdispersion relative to the tree comes from.
width <- runif(n.com, 0.4, 2.5)
regional.abund <- exp(rnorm(n.spp, 1, 0.8))
lam <- outer(seq_len(n.com), seq_len(n.spp), function(i, j) {
  d1 <- (t1[j] - env1[i]) / width[i]
  d3 <- (t3[j] - env2[i]) / 1.5
  regional.abund[j] * exp(-0.5 * d1^2) * exp(-0.25 * d3^2)
})
lam <- lam / rowSums(lam) * runif(n.com, 60, 200)
gradient.sample <- matrix(rpois(length(lam), lam), n.com, n.spp,
                          dimnames = list(rownames(gradient.env),
                                          gradient.phylo$tip.label))
gradient.sample <- as.data.frame(gradient.sample)
stopifnot(all(rowSums(gradient.sample > 0) >= 5), all(colSums(gradient.sample) > 0))

for (o in c("gradient.phylo", "gradient.traits", "gradient.env", "gradient.sample"))
  save(list = o, file = file.path("data", paste0(o, ".rda")), compress = "xz")
summary(rowSums(gradient.sample > 0))
