###############################
# GDS 2023/09/18 Adegenet Tut #
###############################

# Install adegent package
install.packages("adegent", dep=TRUE)

library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

#######################
# CH 3 Object classes #
#######################

# 3.1 genind objects #
data(nancycats)

# logical statemet that checks if this is indeed a genind object
is.genind(nancycats)

nancycats

# this represents the allelic counts (rows 10-18) at the locus (location on the chromosome) "fca8"
# for individual (cat) N224, allele  135 (version of a gene (sequence of nucleotides)) at locus fca8
# is homozygous, meaning that the alleles inherited by the individual at that locus are identical
# whereas for individual N141, allele 129/133 are heterozygous at locus fca8
nancycats[10:18, loc="fca8"]@tab # recommended to use accessors instead of @ or $

# read.genetix reads GENETIX data files and is being put into the object name "obj"
obj <- read.genetix(system.file("files/nancycats.gtx",package = "adegenet"))

# call the R function assigned to the object "obj" using a list to hold the function's arguments instead of writing out the arguments
obj@call 

# evaluates and explicitly states the action performed on the expression, then computes it
toto <- eval(obj@call)

# obj and toto perform the same computation
identical(obj,toto)


# 3.2 genpop objects #
data(nancycats)
catpop <- genind2genpop(nancycats)

catpop

tab(catpop)[1:5, 1:10] # using the accessor for the allele table

# 3.3 Using accessors #
head(indNames(nancycats),10) # first 10 labels of individuals

# append the word cat and the corresponding individual number from nancycats
indNames(nancycats) <- paste("cat",1:nInd(nancycats),sep=".")

# first 10 cat and number labels 
head(indNames(nancycats),10)

locNames(nancycats) # names of loci

# names of loci and allele number 
temp <- locNames(nancycats, withAlleles=TRUE)

# first ten loci and corresponding allele number
head(temp, 10)


obj <- nancycats[sample(1:50,10)]
# accessing the slot 'pop'
pop(obj)

#replicate "newPop" to the first 10 values in the slot "pop"
pop(obj) <- rep("newPop", 10)
pop(obj)

# setting new names for loci
head(colnames(tab(obj)),20)
locNames(obj)

# replace locus name for the first allele
locNames(obj)[1] <- "newLocusName"
locNames(obj)

# now for all 20 loci
head(colnames(tab(obj)),20)


#################################
# CH 4 Importing/Exporting data #
#################################

# 4.1 Importing data frmo GENETIX, STRUCTURE, FSTAT, Genepop #
obj1 <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"))

# you can use this function so that way you don't have to specify the type of file extension
obj2 <- import2genind(system.file("files/nancycats.gtx", package="adegenet"))

all.equal(obj1,obj2)

# 4.2 Importing data from other software #
# .csv files can be converted into genind object using df2genind
# separators can be tricky to use so make sure you know the parameters of the data

# 4.3 Handling presence/absence data #

# 4.4 SNPs data #

# convert SNP to a genind is using df2genind
dat <- matrix(sample(c("a", "t", "g", "c"), 15, replace=TRUE),nrow=3)
rownames(dat) <- paste("genot.", 1:3)
colnames(dat) <- 1:5
dat

obj <- df2genind(dat, ploidy=1)

# obj is a genind containing the SNPs information, which can be used
# for further analysis in adegenet
tab(obj)

################################
# CH 5 Basics of data analysis #
################################

# 5.1 Manipulating the data #
data(microbov)
toto <- genind2genpop(microbov)

toto

popNames(toto)

#subset for the first 3 populations 
titi <- toto[1:3, ]
popNames(titi)

#subset loci using indices or logicals to output locNames
nAll(titi)
tata <- titi[, loc=c(1, 3)]
tata
nAll(tata)

#subset using their explicit name
locNames(titi)
hel5 <- titi[, loc="HEL5"]
hel5
locNames(hel5)

# dropping alleles in the subset
data(nancycats)
nAll(nancycats[1:3, ])
nAll(nancycats[1:3, , drop=TRUE])

# simplify the task of separating data by marker systematically
# using seploc
sepCats <- seploc(nancycats)
class(sepCats)
names(sepCats)
sepCats$fca45

identical(tab(sepCats$fca45), tab(nancycats[,loc="fca45"]))

#separate genotypes in a genind object by population using seppop
data(microbov)
obj <- seppop(microbov)
class(obj)
names(obj)
obj$Borgou
# method to separate data by population and by marker
obj <- lapply(obj,seploc)
names(obj)
class(obj$Borgou)
names(obj$Borgou)
obj$Borgou$INRA63

# method to pool genotypes in different datasets, but having the same marker
# into a single dataset
obj <- seppop(microbov)
names(obj)
newObj <- repool(obj$Borgou, obj$Charolais)
newObj

# the @other slot can be processed during th conversion from genind to genpop
data(sim2pop)
sim2pop
nInd(sim2pop)
head(other(sim2pop)$xy)
dim(other(sim2pop)$xy)
other(genind2genpop(sim2pop, process.other=TRUE))

# 5.2 Using summaries #
toto <- summary(nancycats)
names(toto)
par(mfrow=c(2,2))
plot(toto$n.by.pop, toto$pop.n.all, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")
text(toto$n.by.pop,toto$pop.n.all,lab=names(toto$n.by.pop))
barplot(toto$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")
barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(toto$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)


par(mfrow = c(1,1))
bartlett.test(list(toto$Hexp,toto$Hobs))
t.test(toto$Hexp,toto$Hobs,pair=T,var.equal=TRUE,alter="greater")


# 5.3 testing for hardy-weinberg equilibrium #
library(pegas)
data(nancycats)
cats.hwt <- hw.test(nancycats, B=0)
cats.hwt

# 5.4 Measuring and testing population structure (a.k.a F statistic) #
library("hierfstat")
# This table provides two F statistics F st (pop/total), and F is (ind/pop). These are
# overall measures which take into account all genotypes and all loci.
wc(nancycats)
library(pegas)
ftab <- Fst(as.loci(nancycats))
ftab # per-locus F-stat
colMeans(ftab)

nc <-genind2hierfstat(nancycats)
boot.vc(nc[1], nc[-1])$ci
