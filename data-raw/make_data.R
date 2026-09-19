## Build the package datasets from the raw text files distributed with the
## first edition (Springer supplementary material), kept in inst/extdata so
## that the book can also demonstrate reading them from disk.
library(ape)
ex <- function(f) system.file("extdata", f, package = "lutefisk", mustWork = TRUE)
if (!nzchar(system.file(package = "lutefisk"))) ex <- function(f) file.path("inst/extdata", f)

## The ten-species toy system used in chapters 3 to 8 (and lefse 0.5)
lefse.sample <- read.table(ex("beta.example.sample.txt"), header = TRUE, row.names = 1)
lefse.phylo  <- read.tree(ex("beta.example.phylo.txt"))
lefse.traits <- read.table(ex("beta.example.traits.txt"), header = TRUE, row.names = 1)
stopifnot(identical(read.tree(ex("null.example.phylo.txt"))$tip.label, lefse.phylo$tip.label))
## keep the chapter-specific names as aliases so book code reads naturally
beta.example.sample <- null.example.sample <- partition.example.sample <- lefse.sample
beta.example.phylo  <- null.example.phylo  <- partition.example.phylo  <- comparative.phylo <- lefse.phylo
beta.example.traits <- null.example.traits <- partition.example.traits <- comparative.traits <- lefse.traits

## Chapter 2: a 26-tip tree and a tree with polytomies
my.phylo <- read.tree(ex("my.phylo.newick.file.txt"))
example.poly <- read.tree(ex("example.poly.txt"))

## Chapter 3: PD example (15 species, 5 communities)
PD.example.sample <- read.table(ex("PD.example.sample.txt"), header = TRUE, row.names = 1)
PD.example.phylo  <- read.tree(ex("PD.example.phylo.txt"))

## Chapter 4: FD example (12 species, 5 traits)
FD.example.sample <- read.table(ex("FD.example.sample.txt"), header = TRUE, row.names = 1)
FD.traits <- read.table(ex("FD.traits.txt"), header = TRUE, row.names = 1)

## Chapter 6: 40-tip tree with presence-only communities in three-column form
null.pd.example.phylo <- read.tree(ex("null.pd.example.phylo.txt"))
null.pd.example.sample <- read.table(ex("null.pd.example.sample.txt"), header = FALSE,
	col.names = c("community", "abundance", "species"))

## Chapter 8: environment and coordinates
env.data <- read.table(ex("env.data.txt"), header = TRUE, row.names = 1)
xy.data <- read.table(ex("xy.data.txt"), header = TRUE, row.names = 1)
xy.beta.data <- read.table(ex("xy.beta.data.txt"), header = TRUE, row.names = 1)
new.my.xys <- read.table(ex("new.my.xys.txt"), header = TRUE, row.names = 1)

objs <- c("lefse.sample", "lefse.phylo", "lefse.traits",
	"beta.example.sample", "beta.example.phylo", "beta.example.traits",
	"null.example.sample", "null.example.phylo", "null.example.traits",
	"partition.example.sample", "partition.example.phylo", "partition.example.traits",
	"comparative.phylo", "comparative.traits",
	"my.phylo", "example.poly", "PD.example.sample", "PD.example.phylo",
	"FD.example.sample", "FD.traits", "null.pd.example.phylo", "null.pd.example.sample",
	"env.data", "xy.data", "xy.beta.data", "new.my.xys")
for (o in objs) save(list = o, file = file.path("data", paste0(o, ".rda")), compress = "xz")
