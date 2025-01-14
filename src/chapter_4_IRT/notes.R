# Item Response Theory


# Also known as:
# latent trait analysis
# modern test theory


# Designed for categorical data.
# Characterised by
#* input data is dichotomous or polytomous
#* dimensionality of the underlying trait: uni or multidimensional

# 1. Assessing Dimensionality of a scale
# - categorical principal component analysis
# (Princals, dim. red. method suited for categorical data)
# - exploratory factor analysis (EFA)
# - item factor analysis (IFA)
# (IFA is a variant of EFA for categorical data)

library(MPsychoR)
library(mirt)
library(Gifi)
library(ltm)
library(eRm)
library(mice)

data("RWDQ")
data("zareki")
data("Wilmer")
data("CEAQ")
data("ASTI")
data("WilPat")

zarsub <- zareki[, grep("subtr", colnames(zareki))]

# do items measure a single latent trait or multiple traits?
# let's fit a two-dimensional Princals solution
prinzar <- princals(zarsub)
plot(prinzar, main = "Zareki Loadings")
# multidimensionality = items show in diverging directions

# IFA
# we fit a one and a two factor model and compare the results with a likelihood-ratio (LR) test
fitifa_one <- mirt(zarsub, 1, verbose = FALSE)
fitifa_two <- mirt(zarsub, 2, verbose = FALSE)
anova(fitifa_one, fitifa_two, verbose = FALSE)
plot(fitifa_one)
plot(fitifa_two)

# goals of IRT
# - create scales, gather a set of "high-quality" items
# - score individuals

# Unidimensional, dichotomous IRT models
# Rasch model
# P(X_vi = 1) = (exp(Theta_v + Beta_i)) / (1 + exp(Theta_v + Beta_i))
# we model the p of person "v" scoring 1 on item "i"
# each item gets its item location paramterer Beta_i, which is the item easiness/item difficulty draw_parameter
# each person gets a Theta_v person ability parameter

# the model assumes
# 1., unidimensionality of the trait
# 2., parallel item characteristic curves (ICC)
# 3., local independence
fitrasch1 <- RM(zarsub)
# verifying model fitting by measuring measurement invariance
# this assumes that item paramteres have to be invariant across person subgroups
# best by splitting by binary covariates
timecat <- factor(zareki$time <= median(zareki$time), labels = c("fast", "slow"))
fitLR <- LRtest(fitrasch1, timecat)
# explore which item makes the test fail
Waldtest(fitrasch1, timecat)
plotGOF(fitLR, ctrline = list(col = "gray"), conf = list())
# shows confidence bounds around the diagonal
# 95% confidence ellipses for each item

# test assumptions more explicitly
set.seed(123)
# omitting item 5
T1 <- NPtest(as.matrix(zarsub[, -5]), n = 1000, method = "T1")
T11 <- NPtest(as.matrix(zarsub[, -5]), n = 1000, method = "T11")
# sorted difficulty parameters
fitrasch2 <- RM(zarsub[, -5])
round(sort(-fitrasch2$betapar), 2)

# ICC: item response functions
# the item's behavior along the latent trait
plotjointICC(fitrasch2)
# person paramters
zarppar <- person.parameter(fitrasch2)

# feed back the estimates to the original data and contrast by class
zareki$theta <- zarppar$theta.table[, 1]
summary(aov(theta ~ class, data = zareki))


# one and two-parameter logistic models (1-PL, 2-PL)
# 1-PL estimates a single item discrimination parameter alpha, ICCs parallel
# 2-PL estimates an item discrimination parameter alpha_i for each item, allowing ICCs to cross
# P(X_vi = 1) = (exp(alpha_i(Theta_v - Beta_i))) / (1 + exp(alpha_i(Theta_v - Beta_i)))
# local indenepndence and unidiminsionality are still hold, but the parallel ICCs are relaxed with alpha_i

fit2p12 <- ltm(RWDQ ~ z1) # z1 is a generic placeholder for the single latent dimension
# item parameters
coef(fit2p12) # col1 is the difficulty, col2 is the item discrimination parameters alpha_i, reflecting ICC slopes
# item 22 is weak in difficulty
prinRWDQ <- princals(RWDQ)
plot(prinRWDQ, main = "RWDQ Loadings")

RWDQ <- RWDQ[, -1]
fit2p12 <- ltm(RWDQ ~ z1) # z1 is a generic placeholder for the single latent dimension
# check item fit with Q1
ltm::item.fit(fit2p12)
pdf("pl2.pdf")
plot(fit2p12, items = 1:5, legend = TRUE)
dev.off()
# person parameters
ppars <- ltm::factor.scores(fit2p12, resp.patterns=RWDQ)$score.dat[, "z1"]

# three parameter logistic model (3-PL)
# P(X_vi = 1) = gamma_i + (1 - gamma_i) * (exp(alpha_i(Theta_v - Beta_i))) / (1 + exp(alpha_i(Theta_v - Beta_i)))
# gamma_i reflects probability of 1-response on item "i" due to chance alone, it's called "pseudo-guessing" parameter
# unidimensionality and local indepence are still asssumed

VPMT <- Wilmer[, 3:27]
fit3pl <- tpm(VPMT)
round(head(coef(fit3pl)),3)
pdf("pl3.pdf")
plot(fit3pl, item=1:6, legend=TRUE)
dev.off()


# Unidimensional polytomous IRT models
# Rate Scale Model
# P(X_vi = h) = (exp(h(Theta_v + Beta_i) + omega_h)) / SUM^k ( exp(l(Theta_v + Beta_i) + omega_l))
# theta_v - person parameter
# beta_i - item location parameter (easiness)
# item response categories are h, each gets omega_h (constant across items!) - item differences are reflected by shifts in beta_i
itceaq <- CEAQ[, 1:16] -1
fitrsm <- eRm::RSM(itceaq)
ppar <- eRm::person.parameter(fitrsm)
ifit0 <- eRm::itemfit(ppar)

ind <- match("ceaq10", colnames(itceaq))
itceaq1 <- itceaq[, -ind]
fitrsm1 <- eRm::RSM(itceaq1)
ppar1 <- eRm::person.parameter(fitrsm1)
ifit1 <- eRm::itemfit(ppar1)

ind <- match("ceaq15", colnames(itceaq))
itceaq2 <- itceaq1[, -ind]
fitrsm2 <- eRm::RSM(itceaq2)
ppar2 <- eRm::person.parameter(fitrsm2)
ifit1 <- eRm::itemfit(ppar1)

set.seed(123)
imp <- mice(CEAQ)
gradevec <- complete(imp)$grade
levels(gradevec) <- c("grade56", "grade56", "grade78", "grade78")
eRm::LRtest(fitrsm2, gradevec)
# convert the item parameters into thresholds
thpar <- eRm::thresholds(fitrsm2)

# for polytomous models we get an item-category characteristic curve for each category
pdf("poly1.pdf")
eRm::plotICC(fitrsm2, item=1)
dev.off()
# person-item map
pdf("personitemmap.pdf")
eRm::plotPImap(fitrsm2, latdim="Empathy", main="Person-Item Map CEAQ")
dev.off()

muraki_fit <- mirt::mirt(itceaq1, itemtype="grsm")
plot(muraki_fit)


# partial credit model and generalizaitons
# partial credit model
# we estimate specific item-category parameters
# items don't have to have the same number of categories
# P(X_vih = 1) = (exp(h*Theta_v + Beta_ih)) / SUM^k ( exp(l * Theta_v + Beta_ih))
# it differs from RSM, that each item-category gets is own item-category paramter Beta_ih
PGitems <- ASTI[, c(11, 14, 15, 17, 18, 23)]
fitPCM <- eRm::PCM(PGitems)
threPCM <- eRm::thresholds(fitPCM)
eRm::plotPImap(fitPCM, latdim="Presence/Growth", main="Person-Item Map ASTI")
eRm::plotICC(fitPCM, item=2)

#generalized partial credit model
# Muraki generalizes the model by adding an item discrimination parameter alpa_i
# P(X_vih = 1) = (exp(alpha_i(h*Theta_v + Beta_ih))) / SUM^k ( exp(alpha_i(l * Theta_v + Beta_ih)))
# slopes vary across items, contstant within item
STitems <- ASTI[, c(2, 4, 7, 13, 16, 24, 25)]
stpcm <- ltm::gpcm(STitems, constraint = "rasch") #pcm
stgpcm <- ltm::gpcm(STitems) # gpcm
anova(stpcm, stgpcm)

# Graded Response Model
# P(X_vih >= 1) = (exp(alpha_i(Theta_v - Beta_ih))) / (1 +  exp(alpha_i(Theta_v - Beta_ih)))
# formulates P for category h or higher
# beta_ih is called the category boundary location
# items are allowed to have different number of categories
fitgrm <- ltm::grm(STitems)
ppargrm <- ltm::factor.scores(fitgrm)

#operation characteristic curve
pdf("grm1.pdf")
plot(fitgrm, items=1, type="OCCu")
dev.off()
pdf("grmICC.pdf")
plot(fitgrm, items=1, type="ICC")
dev.off()

# Nominal Response Model
# this model abandons the requirement of ordinality in item categories
# ideal for multiple choice or unordered responses
# each item-category gets its own individual discrimination parameter beside the item-category location
# P(X_vih = h) = (exp(alpha_ih(Theta_v - Beta_ih))) / (1 +  exp(alpha_ih(Theta_v - Beta_ih)))
# alpha_ih is the item-category discrimination parameter
# ICCs don't have to be parallel across or within items
wpit15 <- WilPat[,1:15]
wpiprin <- Gifi::princals(wpit15,  ordinal=FALSE)
plot(wpiprin)
elim <- c("Nationalism", "Patriotism", "ChurchAuthority", "Obedience")
ind <- match(elim, colnames(wpit15))
wpitnew <- wpit15[, -ind]
wpihom <- Gifi::homals(wpitnew)
plot(wpihom)

fitnrm <- mirt::mirt(wpitnew, 1, itemtype="nominal")
mirt::itemplot(fitnrm, item=5)
#goodness-of-fit
mirt::M2(fitnrm)

# Item and Test Information
# in which area of the trait an item is particularly informative
# in what is the defree to which an item reduces the uncertainty in estimation of a person's trait value
# item information
plot(fitnrm, type="infotrace", main="Item Information")

# rumination data
save_plot <- function(p, plotname){ 
   pdf(plotname)
   p
   dev.off()
}
rumination  <- readRDS("~/repos/rumination_lasso/data/data.RDS")
scrs_ind <- match(paste0("SCRS_",1:10), colnames(rumination))
scrs  <- rumination[,scrs_ind]
scrs_fit <- mirt::mirt(scrs, 1)
pdf("SCRSfit.pdf")
plot(scrs_fit, type="infotrace", main="Item Information Curves for SCRS")
dev.off()
plot(scrs_fit, type="info", main="Self Critical Rumination Scale Test Information")

rrs_ind <- match(paste0("RRS_", 1:10), colnames(rumination))
rrs <- rumination[,rrs_ind]
rrs_fit <- mirt::mirt(rrs, 1)
rrsplot <- plot(rrs_fit, type="infotrace", main="Item Information Curvers for RRS")
save_plot(plot(rrs_fit, type="infotrace"), "RRSfit.pdf")
plot(rrs_fit, type="info", main="Ruminative Response Scale Test Information")


